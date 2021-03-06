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

        xmer_seq_dict = {}

        # create list of Xmer sequences
        for index in range( len( sequences ) ):

            name, sequence = oligo.subset_lists_iter( names[ index ], sequences[ index ], options.XmerWindowSize, options.stepSize )

            for index in range( len( sequence ) ):
                if oligo.is_valid_sequence( sequence[ index ], options.minLength, options.percentValid ):
                    value = [ options.redundancy, name[ index ] ]
                    xmer_seq_dict[ sequence[ index ] ] = value

        # create dict of Ymer sequences
        ymer_seq_dict = {} 

        # Break each ymer up into subsets of xmer size
        for index in range( len( sequences ) ):

            name, sequence = oligo.subset_lists_iter( names[ index ], sequences[ index ], options.YmerWindowSize, options.stepSize )


            for index in range( len( sequence ) ):

                if oligo.is_valid_sequence( sequence[ index ], options.minLength, options.percentValid ):
                    ymer_seq_dict[ sequence[ index ] ] = name[ index ]

        total_ymers = len(ymer_seq_dict)

        array_design = {}
        array_xmers = {}
        to_add = []
        ymer_xmers = []
        iter_count = 0

        while True:
            #reset max score at the beginning of each iteration
            max_score = 0
            for current_ymer in ymer_seq_dict.keys(): 
                # calculate the score of this ymer
                score, subset_ymer = calculate_score( current_ymer, xmer_seq_dict, options.XmerWindowSize, 1 )
            
                if score > max_score:
                    to_add = list()
                    max_score = score
                    to_add.append( current_ymer )
                    ymer_xmers = [subset_ymer]
                elif score == max_score:
                    to_add.append( current_ymer )
                    ymer_xmers.append(subset_ymer)

            random_index = random.choice( range( len( to_add ) ) )
            oligo_to_remove = to_add[ random_index ]
            chosen_xmers = ymer_xmers[ random_index ]
    #        array_xmers.update(chosen_xmers)
            for each in chosen_xmers:
                array_xmers[each] = array_xmers.get(each, 0)+1

            # subtract from the score of each xmer within the chosen ymer
            for item in chosen_xmers:
                if item in xmer_seq_dict:
                # We dont' want negative scores
                    if xmer_seq_dict[ item ][ 0 ] > 0:
                        xmer_seq_dict[ item ][ 0 ] -= 1
                else: print ( "%s - not found in xmer dict!!!" % (item) )

            iter_count += 1

            if len( ymer_seq_dict ) == 0 or max_score <= 0:
                print ( "Final design includes %d %d-mers (%.1f%% of total) " % (len(array_design), options.YmerWindowSize, (len(array_design)/float(total_ymers))*100) )
    #            average_redundancy = sum( xmer_seq_dict[ item ][ 0 ] for item in xmer_seq_dict ) / len( xmer_seq_dict )
                print( "%d unique %d-mers in final %d-mers (%.2f%% of total)" % (len(array_xmers), options.XmerWindowSize, options.YmerWindowSize, (float(len(array_xmers))/len(xmer_seq_dict))*100))
                print( "Average redundancy of %d-mers in %d-mers: %.2f" % (options.XmerWindowSize, options.YmerWindowSize, sum(array_xmers.values())/float(len(array_xmers))))
                if len(array_design)<min_ymers:
                    min_ymers = len(array_design)
                    best_xmer_seq_dict = xmer_seq_dict
                    del(xmer_seq_dict)
                    best_array_design = array_design
                    del(array_design)
                break
        

            try:
                array_design[ oligo_to_remove ] = ymer_seq_dict[ oligo_to_remove ]
                del ymer_seq_dict[ oligo_to_remove ]
            except KeyError:
                continue

            if not iter_count % 250:
                print( "Current Iteration: " + str( iter_count ) ) 
    #            print( "Number of output ymers: " + str( len( array_design ) ) )
                print( "Current xmer dictionary score: " + str( sum( item[ 0 ] for item in xmer_seq_dict.values() ) ) )



    write_outputs( best_xmer_seq_dict, options.outPut )

    names = []
    sequences = []

    # Write resulting oligos to file
    for sequence, name in best_array_design.items():
        names.append( name )
        sequences.append( sequence )

    oligo.write_fastas( names, sequences, output_name=options.outPut + "_R" + str( options.redundancy ) + ".fasta" ) 


def calculate_score( ymer, comparison_dict, window_size, step_size ):
    """
        Calculates the score of a ymer
    """
    name, subset_ymer = oligo.subset_lists_iter( "", ymer, window_size, step_size )
    total = 0
    for current_ymer in subset_ymer:
       if current_ymer in comparison_dict and isinstance( comparison_dict[ current_ymer ], list ):
           total +=  comparison_dict[ current_ymer ][ 0 ] 
    return total, subset_ymer
    
def write_outputs( seq_dict, out_file ):
    """
        Writes tab delimited key-value pairs to specified output file
        Each line consists of key \t value  
    """
    file = open( out_file, 'w+' )
    for xmer, score in seq_dict.items():
        file.write( xmer + '\t' + str( score[ 0 ] ) + '\n' )
    file.close()
                


    
     


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
