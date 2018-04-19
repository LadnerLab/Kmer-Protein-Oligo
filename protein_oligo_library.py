#!/usr/bin/env python3
   
def read_fasta_lists( file_to_read ):
    """
       Reads a list of fastas from file_to_read
    
       Returns:
        names- a list of names of the sequences found in the fasta file
        sequences- a list of the sequences found in the fasta file
    """

    file_in = open( file_to_read, 'r' )
    count = 0

    names = []
    sequences = []
    current_sequence = ''

    for line in file_in:
        line = line.strip()
        if line and line[ 0 ] == '>':
            count += 1
            names.append( line[ 1: ] )
            if count > 1:
                sequences.append( current_sequence )
            current_sequence = ''

        else:
            current_sequence += line

    sequences.append( current_sequence )
    file_in.close()

    return names, sequences
 
def write_fastas( names_list, sequence_list, output_name="out.txt" ):
    """
        Writes a fasta file from a list of names and sequences to output file provided

    """
    out_file = open( output_name, 'w+' )
    for index in range( len( names_list ) ):
        out_file.write( '>' + names_list[ index ] + '\n' + 
                        sequence_list[ index ] + '\n'
                      )
    out_file.close()
        

def char_in_string( test_string, character ):
    """
        Checks if a character is found within a given string
    """
    for index in range( len( test_string ) ):
        if test_string[ index ] == character:
            return True
    return False

def percentage_of_char_in_string( test_string, character ):
    """
        Calculates what percent of a given test_string is given character
    
        Params:
           test_string- string to test for character
           character- character to test for in string
       Returns:
           floating point value percent of the string that
           is character
    """
    length = len( test_string )
    char_count = 0.0

    for index in range( length ):
        if test_string[ index ] == character:
            char_count += 1
    return ( char_count / length ) * 100

def count_char_in_string( test_string, character ):
    """
        Counts how much of a certain character are in test_string

        Params:
            test_string- string to count
            character- character to count in string
        Returns:
            integer value representing the number of character 
                were found in test_string
    """
    length = len( test_string )
    count = 0

    for index in range( length ):
        if test_string[ index ] == character:
            count += 1
    return count

def min_concurrent_chars( test_string, delimeter_char ):
    """
        Finds the minimum number of concurrent non-delimiter char in 
        test_string
    
        Params:
          test_string- string to search
          delimeter_char- character to reset the count
        Returns:
          integer value, the smallest amount of concurrent characters
               between delimeter character
    """

    split_string = test_string.split( delimeter_char )
    min_length = len( split_string[ 0 ] )

    for substring in split_string[ 1: ]:
        current_length = len( substring )
        if current_length > 0 and current_length < min_length:
            min_length = current_length

    return min_length


def remove_char_from_string( test_string, to_remove ):
    """
        Removes character to_remove from string test_string
        Note: Case sensitive method, 'a' != 'A'
        Returns:
           test_string, minus any instance of to_remove character
    """
    output = ""
    for index in range( len( test_string ) ):
        if test_string[ index ] != to_remove:
            output += test_string[ index ]
    return output


def create_list_of_uniques( names, sequences ):
    """
       Removes duplicates from a list
       Params:
          names- names of the sequences
       sequences:
          a list of sequences who may or may not be unique
    """
    return_names = []
    return_sequences = []

    unique_values = set()

    for index in range( len( sequences ) ):
        starting_length = len( unique_values )
        unique_values.add( sequences[ index ] )

        if len( unique_values ) > starting_length:
            return_names.append( names[ index ] )
            return_sequences.append( sequences[ index ] )
    return return_names, return_sequences

def create_valid_sequence_list( names_list, sequence_list, min_length, percent_valid, start, end ):
   """
       Creates a sequence list of valid sequences.
       A valid sequence is defined by not having any 'X' characters,
       and not violating the parameters of either min_length or percent_valid 
       
       Returns:
           a list of names of those sequences that were valid, with the new bounds appended to the name
           a list of the sequences that were found valid
   """
   valid_names = []
   valid_sequences = []

   for sequence in range( len( sequence_list ) ):
      current_sequence = sequence_list[ sequence ][ start : end ]

      if is_valid_sequence( current_sequence, min_length, percent_valid ):
           valid_names.append( names_list[ sequence ] )
           current_sequence = remove_char_from_string( current_sequence, '-' )
           valid_sequences.append( current_sequence[ start : end ] )

   names_list = append_suffix( valid_names, start, end )

   return names_list, valid_sequences

def is_valid_sequence( sequence, min_length, percent_valid ):
   """
       Determines whether a given sequence is valid 
       A valid sequence is defined by not having any 'X' characters,
           and not violating the parameters of either min_length or percent_valid 
   """
   if not char_in_string( sequence, 'X' ):
       if min_length is None:
           return percentage_of_char_in_string( sequence, '-' ) < ( 100 - percent_valid )
       else:
           return ( min_concurrent_chars( sequence, '-' ) >= min_length )
   return False

         
def append_suffix( string, start, end ):
   """
       Appends _start_end to a string
   """
   return "%s_%s_%s" % ( string, str( start ), str( end ) ) 


def subset_lists_iter( name, sequence, window_size, step_size ):
    new_names = []
    new_seqs = []

    start = 0
    end = window_size

    while end <= len( sequence ):
       new_seqs.append( sequence[ start : end ] )
       new_names.append( append_suffix( name, start + 1, end ) )

       start += step_size
       end = start + step_size + window_size - 1

    return new_names, new_seqs

def subset_lists( name, sequence, window_size, step_size ):
   """
       Creates a list of subsets of windowSize size in intervals of stepSize
       Note: Uses recursive subset_lists_helper for operations
   
       Params:
            name: String name of sequence to be split up
            sequence: String sequence to be split up into a list
       Returns:
            a list of the split up names, with a suffix applied, and a list of the segments 
            of the list, as specified
   """
   new_names = []
   new_seqs = []
   return subset_lists_helper( name, sequence, new_names, new_seqs, window_size, step_size, 0, window_size )

def subset_lists_helper( name, sequence, name_arr, seq_arr, window_size, step_size, start, end ):
    """
        Recursive helper method called by subset_lists
    """
    if start + window_size < len( sequence ):
       if len( sequence[ start: end ] ) == 1:
           return
       seq_arr.append( sequence[ start : end ] ) 
       name_arr.append( append_suffix( name, start + 1, end ) )

       subset_lists_helper( name, sequence, name_arr, seq_arr, window_size, step_size, start + step_size, start + step_size + window_size )
    return name_arr, seq_arr
   


            
    
