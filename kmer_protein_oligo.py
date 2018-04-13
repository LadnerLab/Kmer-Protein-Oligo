#! /usr/bin/env python3
import optparse
import sys
import protein_oligo_library as oligo


def main():
    usage = "usage: %prog [options]"
    option_parser = optparse.OptionParser( usage )
    
    add_program_options( option_parser )
    
    step_size = 1


def add_program_options( option_parser ):
   option_parser.add_option( '-q', '--query', help = "Fasta query file of sequence to be used by program. [None, Required]"
   )

   option_parser.add_option( '-w', '--windowSize', type = 'int', \
                             default = 100, \
                             help = "Amount of characters from each alignment sequence to look at. [100]"
   )
   option_parser.add_option( '-o', '--outPut', default = "oligo_out.txt", help = "Name of file program output will be written to. [oligo_out.txt]"
   )

   option_parser.add_option( '-r', '--redundancy', type = 'int', default = 8, help = "A number specifying the redundancy to be used to each kmer [8]" )



if __name__ == '__main__':
    main()
