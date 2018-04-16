#! /usr/bin/env python3
import optparse
import sys
import protein_oligo_library as oligo


def main():
    usage = "usage: %prog [options]"
    option_parser = optparse.OptionParser( usage )
    
    add_program_options( option_parser )

    options, arguments = option_parser.parse_args()

    names, sequences = oligo.read_fasta_lists( options.query )
    xmer_seq_dict = {}

    # create list of Xmer sequences
    for index in range( len( sequences ) ):

        # We know the recursion is finite becuase the sequence is finite
        sys.setrecursionlimit( len( sequences[ index ] ) + 50 )
        name, sequence = oligo.subset_lists( names[ index ], sequences[ index ], options.XmerWindowSize, options.stepSize )

        for sub_sequence in sequence:
            if oligo.is_valid_sequence( sub_sequence, options.minLength, options.percentValid ):
                xmer_seq_dict[ sub_sequence ] = options.redundancy


    ymer_seq_dict = {}

    # create list of Ymer sequences
    for index in range( len( sequences ) ):

        # We know the recursion is finite becuase the sequence is finite
        sys.setrecursionlimit( len( sequences[ index ] ) + 50 )
        name, sequence = oligo.subset_lists( names[ index ], sequences[ index ], options.YmerWindowSize, options.stepSize )

        for sub_sequence in sequence:
            if oligo.is_valid_sequence( sub_sequence, options.minLength, options.percentValid ):
                ymer_seq_dict[ sub_sequence ] = options.redundancy
    print( ymer_seq_dict )




    
     


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

   option_parser.add_option( '-r', '--redundancy', type = 'int', default = 8, help = "A number specifying the redundancy to be used to each kmer [3]" )

   option_parser.add_option( '--stepSize', type = 'int', default = 1, help = (
      "Step size to move over after each subset of windowSize characters has been read"
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
