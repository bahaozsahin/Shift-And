from Bio import SeqIO
import numpy as np
from numpy.core.fromnumeric import shape
from scipy.ndimage.interpolation import shift
import time

#indirilmiş fn uzantılı dosyanın açılması
input_file = open('gene.fna', 'r')


#Eşleştirilecek sekanslar tanımlandı
seq_input_first = 'GAATTCATTTCTCCCGCTGCCCCATCTCT'
seq_input_second = 'TCACG'
seq_input_third = 'TCTCG'

#Eşleşme yapılacak arraylerin U(A) U(T) U(G) U(C) vektörlerinin tanımlanması

def generate_u(sequence_input):
    vector_a = np.zeros(len(sequence_input), dtype= int)
    vector_t = np.zeros(len(sequence_input), dtype= int)
    vector_g = np.zeros(len(sequence_input), dtype= int)
    vector_c = np.zeros(len(sequence_input), dtype= int)
    index = 0
    for letter in sequence_input:
        if(letter == 'A'):
            vector_a[index] = 1
        if(letter == 'T'):
            vector_t[index] = 1
        if(letter == 'G'):
            vector_g[index] = 1
        if(letter == 'C'):
            vector_c[index] = 1
        index += 1

    return (vector_a,vector_t, vector_g, vector_c)


#Bitshift operasyonu için array bir adım sağa kaydırıldıktan sonra en başına 1 eklendi
def bitshift(vector_input):
    result = shift(vector_input, 1, cval=1)
    return result


def check_matches_using_and(bitshifted_vector_of_prev, current_u_of_column):
    result = np.zeros(len(bitshifted_vector_of_prev), dtype= int) 
    for i in range(0, len(bitshifted_vector_of_prev)):
        result[i] = bitshifted_vector_of_prev[i]&current_u_of_column[i]
    return result


def checking_matches_in_matrix(mask_matrix):
    p,t = np.shape(mask_matrix)
    counter_of_matches = 0
    #print(p, t)

    for i in range(0, t):
        if(mask_matrix[0][i] == 1):
            counter_for_matches = 0

            for j in range(0, p):
                #print(i, j)
                if((i+j) < t):
                    #print(i+j)
                    if(mask_matrix[j][j+i] == 1):
                        counter_for_matches += 1
                    if(counter_for_matches == p):
                        #print(j, j+i)
                        counter_of_matches += 1
                        #print(counter_of_matches)
    if(counter_of_matches == 0):
        print('Bulunamadi')
    print(counter_of_matches)
    return counter_of_matches


def main():
    for sequence in SeqIO.parse('gene.fna','fasta'):
        seq_t = sequence.seq
    
    vector_a, vector_t, vector_g, vector_c = generate_u(seq_input_first)
    mask_matrix = np.zeros(shape = (len(seq_input_first), len(seq_t)), dtype= int)
    index = 0
    
    initializing_column = np.zeros(len(seq_input_first), dtype= int)
    prev_column = initializing_column
    bitshifted_prev_column = bitshift(prev_column)

    for letter in seq_t:    
        if(letter == 'A'):
            mask_vector = check_matches_using_and(bitshifted_prev_column, vector_a)   
            for i in range(0,len(mask_vector)):
                mask_matrix[i, index] = mask_vector[i]
            prev_column = mask_vector
        if(letter == 'T'):
            mask_vector = check_matches_using_and(bitshifted_prev_column, vector_t)   
            for i in range(0,len(mask_vector)):
                mask_matrix[i, index] = mask_vector[i]
            prev_column = mask_vector
            
        if(letter == 'G'):
            mask_vector = check_matches_using_and(bitshifted_prev_column, vector_g)   
            for i in range(0,len(mask_vector)):
                mask_matrix[i, index] = mask_vector[i]
            prev_column = mask_vector
            
        if(letter == 'C'):
            mask_vector = check_matches_using_and(bitshifted_prev_column, vector_c)   
            for i in range(0,len(mask_vector)):
                mask_matrix[i, index] = mask_vector[i]
            prev_column = mask_vector
        index += 1
        bitshifted_prev_column = bitshift(prev_column)
    
    
    #print(mask_matrix)
    checking_matches_in_matrix(mask_matrix)

if __name__ == '__main__':
    start = time.time()
    main()
    end = time.time()
    difference= end - start
    print("The time between start and end of the function is (in seconds)= ", difference)
