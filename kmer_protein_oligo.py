#! /usr/bin/env python3
import optparse
import sys
import protein_oligo_library as oligo
import random


def main():
    usage = "usage: %prog [options]"
    option_parser = optparse.OptionParser( usage )
    
    add_program_options( option_parser )

    options, arguments = option_parser.parse_args()

    names, sequences = oligo.read_fasta_lists( options.query )

    min_ymers = 999999999999999999999999999999999

    for i in range(options.iterations): 

        #Values are now all of the locations of the xmer in the input seqs
        xmer_seq_dict = {}

        # create list of Xmer sequences
        for i in range( len( sequences ) ):
            xdict = oligo.subset_lists_iter( i, sequences[i], options.XmerWindowSize, options.stepSize )
            for x, locs in xdict.items():
                #????? Is it appropraite to use the same minLength and percentValid for xmers and ymers????
                if oligo.is_valid_sequence( x, options.minLength, options.percentValid ):
                    if x in xmer_seq_dict: 
                        xmer_seq_dict[x].update(locs)
                    else: xmer_seq_dict[x] = set(locs)

        print( "xmer_dict done!" )

        # create dict of Ymer sequences
        ymer_name_dict = {} 
        ymer_loc_dict = {} 

        # Break each sequence up into ymers
        for i in range(len(sequences)):
            ydict = oligo.subset_lists_iter( names[i], sequences[i], options.YmerWindowSize, options.stepSize )

            #For each ymer, break into xmers and record locations of those ymers
            for y, locs in ydict.items():
                if y not in ymer_loc_dict and oligo.is_valid_sequence( y, options.minLength, options.percentValid ):
                    ymer_name_dict[y] = locs[0]
                    ymer_loc_dict[y] = oligo.component_xmer_locs(y, xmer_seq_dict, options.XmerWindowSize, options.stepSize)

        print( "ymer_dict done!" )

        array_design = {}
        array_xmers = {}
        to_add = []
        iter_count = 0

        while True:
            #reset max score at the beginning of each iteration
            max_score = 0
            
            for current_ymer, locs in ymer_loc_dict.items(): 
                # calculate the score of this ymer
                if len(locs) > max_score:
                    max_score = len(locs)
                    to_add=[current_ymer]
                elif len(locs) == max_score:
                    to_add.append( current_ymer )

            #Terminate the while loop if any of these conditions are met
            if len( ymer_loc_dict ) == 0 or max_score <= 0 or len(array_xmers)/float(len(xmer_seq_dict))>=options.minXmerCov:
                total_ymers = len(ymer_loc_dict)
                print ( "Final design includes %d %d-mers (%.1f%% of total) " % (len(array_design), options.YmerWindowSize, (len(array_design)/float(total_ymers))*100) )
                print( "%d unique %d-mers in final %d-mers (%.2f%% of total)" % (len(array_xmers), options.XmerWindowSize, options.YmerWindowSize, (float(len(array_xmers))/len(xmer_seq_dict))*100))
                average_redundancy = sum(array_xmers.values()) / len(xmer_seq_dict)
                print( "Average redundancy of %d-mers in %d-mers: %.2f" % (options.XmerWindowSize, options.YmerWindowSize, sum(array_xmers.values())/float(len(array_xmers))))

                #To log info about the best iteration
                if len(array_design)<min_ymers:
                    min_ymers = len(array_design)
                    best_xmer_seq_dict = xmer_seq_dict
                    del(xmer_seq_dict)
                    best_array_design = array_design
                    del(array_design)
                break

            #Randomly choose a ymer from the list of ymers with equal coverage, get info and remove from ymer dict
            oligo_to_remove = to_add[random.choice(range(len(to_add)))]
            covered_locs = ymer_loc_dict[oligo_to_remove]
            del(ymer_loc_dict[oligo_to_remove])
            
            #Add info about the contained xmers to the array_xmer dict
            these_xmers = oligo.subset_lists_iter('', oligo_to_remove, options.XmerWindowSize, 1)
            for x in these_xmers:
                array_xmers[x] = array_xmers.get(x,0)+1
            
            #Remove covered xmer locations from the remaining ymers
            for k in ymer_loc_dict.keys():
                ymer_loc_dict[k] = ymer_loc_dict[k].difference(covered_locs)
            
            iter_count += 1
            array_design[oligo_to_remove] = ymer_name_dict[ oligo_to_remove ]

            if not iter_count % 100:
                print( "Current Iteration: " + str( iter_count ) ) 
#                print( "Number of output ymers: %d" % len(array_design))
                print( "Current max score: %d" % (max_score))

#    write_outputs( best_xmer_seq_dict, options.outPut )

    names = []
    sequences = []

    # Write resulting oligos to file
    for sequence, name in best_array_design.items():
        names.append( name )
        sequences.append( sequence )

    oligo.write_fastas( names, sequences, output_name=options.outPut + "_R" + str( options.redundancy ) + ".fasta" ) 


# def calculate_score( ymer, comparison_dict, window_size, step_size ):
#     """
#         Calculates the score of a ymer
#     """
#     name, subset_ymer = oligo.subset_lists_iter( "", ymer, window_size, step_size )
#     total = 0
#     for current_ymer in subset_ymer:
#        if current_ymer in comparison_dict and isinstance( comparison_dict[ current_ymer ], list ):
#            total +=  comparison_dict[ current_ymer ][ 0 ] 
#     return total, subset_ymer
#     
# def write_outputs( seq_dict, out_file ):
#     """
#         Writes tab delimited key-value pairs to specified output file
#         Each line consists of key \t value  
#     """
#     file = open( out_file, 'w+' )
#     for xmer, score in seq_dict.items():
#         file.write( xmer + '\t' + str( score[ 0 ] ) + '\n' )
#     file.close()
                


    
     


def add_program_options( option_parser ):
   option_parser.add_option( '-q', '--query', help = "Fasta query file of sequence to be used by program. [None, Required]"
   )

   option_parser.add_option( '-x', '--XmerWindowSize', type = 'int', \
                             default = 100, \
                             help = "Amount of characters from each Xmer alignment sequence to look at. [100]"
   )

   option_parser.add_option( '-y', '--YmerWindowSize', type = 'int', \
                             default = 100, \
                             help = "Amount of characters from each Ymer alignment sequence to look at. [100]"
   )

   option_parser.add_option( '-o', '--outPut', default = "oligo_out.txt", help = "Name of file program output will be written to. [oligo_out.txt]"
   )

   option_parser.add_option( '-r', '--redundancy', type = 'int', default = 3, help = "A number specifying the redundancy to be used to each kmer [3]" )

   option_parser.add_option( '-i', '--iterations', type = 'int', default = 1, help = "Number of independent iterations to run. The result with the fewest oligos will be output [1]" )


   option_parser.add_option( '--stepSize', type = 'int', default = 1, help = (
      "Step size to move over after each subset of windowSize characters has been read [1]"
      )
      )

   option_parser.add_option( '-c', '--minXmerCov', type = 'float', default = 1, help = "Minimum proportion of total xmers that need to be covered in design for process to terminate [1]" )

   option_parser.add_option( '-l', '--minLength', type = 'int', help = (
      "Minimum length of concurrent non-dash characters that must be present in order for "
      "the sequence to be considered valid, sequences with a maximum length of concurrent non-dash "
      "characters less than this parameter will not be included in program output. [None, Required] "
   )
   )
   option_parser.add_option( '-p', '--percentValid', type = 'float', default = 90.0, help = (
      "Percent of non '-' characters present in order for the sequence to be considered valid, "  
      "sequences with less than specified amount will not be present in program out put. [90.00] "
   )
   )
 
 



if __name__ == '__main__':
    main()
