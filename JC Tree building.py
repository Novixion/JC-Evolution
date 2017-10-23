print('Welcome to JC Evolution Tree Generation') # this is the title

import random #idk what these do but i need a random number generator
import numpy
import scipy
import math

randomtree = input('Would you like to generate a random tree? (yes/no) ')
#notes are written using pounds
#int() is for integer input"

lengt = int(input('How long would you like your sequence to be? '))
term = int(input('How many sequences do you want back? '))
alpha = float(input('Please input the mutations per site per generation: '))
#total terms are term for end, lengt for length of sequence and randomtree for gen rand or not

#generate or input the initial sequence and store it in isq
if randomtree == "yes" or randomtree == "Yes":
    isq = ''
    for m in range(lengt):
        letterm = random.choice(['A','T','G','C'])
        isq += letterm
        #this stores the sequence in isq
        #generates the DNA sequence
    print('Generating random tree')
elif randomtree == "no" or randomtree == "No":
    print('Generating tree with given sequence')
    isq = input('Please input the DNA sequence ') #input sequence           
else:
    print('CODE ERROR BEEP BOOP BEEP BOOP FORCE QUIT PLS')
    print('Please try again using only yes or no')

n = 1
time = 0
D = numpy.zeros((term+1,term+1))

print('This is the sequence you are using:', isq)
           
#doing the JC evolution

Mi = numpy.matrix([[1 - alpha, alpha / 3, alpha / 3, alpha / 3],
                  [alpha / 3, 1 - alpha, alpha / 3, alpha / 3],
                  [alpha / 3, alpha / 3, 1 - alpha, alpha / 3],
                  [alpha / 3, alpha / 3, alpha / 3, 1 - alpha]])
        #INSERT MATRIX HERE
        #JC Matrix
M = Mi
slist = [isq]
timelist = [0]

while n <= term:
    #read isq
    time = time + 1
    nsq = ''
    for letter in isq:
        ranmut = random.uniform(0,1)        #generate a random number between 0 and 1

        if letter == 'A':
            if ranmut <= M[0,0]:        #first term in the matrix
                letter = 'A'            #prob is A
            #somehow replace A with A in the original sequence
            elif ranmut <= M[0,1] + M[0,0]:
                letter = 'T'
            #replace A with T
            elif ranmut <= M[0,2] + M[0,1] + M[0,0]:
                letter = 'G'
            #replace A with G
            elif ranmut <= 1:
                letter = 'C'
            #replace A with C
        
        
        elif letter == 'T':         
            if ranmut <= M[1,0]:
                letter = 'A'
            #somehow replace c with A in the original sequence
            elif ranmut <= M[1,1] + M[1,0]:
                letter = 'T'
            #replace T with T
            elif ranmut <= M[1,2] + M[1,1] + M[1,0]:
                letter = 'G'
            #replace T with G
            elif ranmut <= 1:
                letter = 'C'
            #replace T with C


        elif letter == 'G':
            if ranmut <= M[2,0]:
                letter = 'A'
            #somehow replace G with A in the original sequence
            elif ranmut <= M[2,1] + M[2,0]:
                letter = 'T'
            #replace G with T
            elif ranmut <= M[2,2] + M[2,1] + M[2,0]:
                letter = 'G'
            #replace G with G
            elif ranmut <= 1:
                letter = 'C'
            #replace G with C


        elif letter == 'C':
            if ranmut <= M[3,0]:
                letter = 'A'
            #somehow replace C with A in the original sequence
            elif ranmut <= M[3,1] + M[3,0]:
                letter = 'T'
            #replace C with T
            elif ranmut <= M[3,2] + M[3,1] + M[3,0]:
                letter = 'G'
            #replace C with G
            elif ranmut <= 1:
                letter = 'C'
            #replace C with C
        nsq += letter

        
    M = M.dot(Mi)    #move the matrix one more dimension down
    
#now check for if the sequence is the same if not add to the list, if it is then ignore and repeat
    if isq != nsq:
        slist.append(nsq)
        isq = nsq
        timelist.append(time)
        n = n + 1
        M = Mi

#create the Jukes Cantor Distance Matrix
dis = D
dif = 0
for l in range(0,term+1):
    for k in range(0,term+1):
        dif = 0
        for letter, let in zip(slist[l], slist[k]):
            if letter != let:
                dif = dif + 1
        dis[l,k] = dis[k,l] = dif/lengt

print('Proportion different = ', dis)

for i in range(1,term+1):
    for j in range(1,term+1):
        D[i,j] = D[j,i] =  -3/4*numpy.log(1-4/3*dis[i,j])
        
#final step list the sequences and the final distances
print('')
num = 1
for sequence in slist:
    print(num, sequence)
    num = num + 1
print('')
print('Mutation times: ', timelist)
print('')
print('JC Distance Matrix in Order 1 -', term)
print(D)



















#Shoulda just done this in MatLab
