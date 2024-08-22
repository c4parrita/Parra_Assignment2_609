#!/usr/bin/env python
# coding: utf-8

# In[10]:


def SuffixArray (genomef, readf):
    
    gfile = open(genomef, "r")
    next(gfile)
    
    #concatenate lines into 1 genome string
    genome=''
    for l in gfile:
        l = l.strip()
        genome = genome+l
        
    #create suffix array with all indicies     
    sufarr = list(range(len(genome)))
    
    #organize suffix array lexographically 
    sufarr.sort(key=lambda x: genome[x:])

    #obtain read from file 
    rfile = open(readf, "r")
    next(rfile)
    read = next(rfile).strip()
    
    occur = []
    
    #binary seacrch for read (bisection)
    low = 0
    high = len(sufarr) - 1
    
    while low <= high:
        
        #find new midpoint based on low and high ranges 
        midpoint = (high+low) // 2
        
        #if occurence has been found at midpoint add to array
        if genome[sufarr[midpoint]:sufarr[midpoint]+len(read)] == read:
            occur.append(sufarr[midpoint])
            
            #look at the right side for more occurences
            i = midpoint + 1
            while i < len(sufarr):
                if genome[sufarr[i]:sufarr[i]+len(read)] == read:
                    occur.append(sufarr[i])
                    i+=1
                    
                #if suffix does not match pattern, no more occurences to the right 
                else:
                    break
                    
            #look at the left side for more occurences
            j = midpoint - 1
            while j >= 0:
                if genome[sufarr[j]:sufarr[j]+len(read)] == read:
                    occur.append(sufarr[j])
                    j-=1
                    
                #if suffix does not match pattern, no more occurences to the left 
                else:
                    break
            occur = sorted(occur)        
            print(occur)
            with open('occurences.txt', 'w') as file:
  
                for oc in occur:
                    file.write(str(oc) +'\n')
                    
            return 1
        
        #occurence most likely to the left of current midpoint  
        elif genome[sufarr[midpoint]:sufarr[midpoint]+len(read)] < read:
            low = midpoint + 1
        
        #occurence most likely to the right of current midpoint  
        else:
            high = midpoint - 1
  
    return 1


# In[11]:


SuffixArray('Assignment2_refgenome.fasta', 'Assignment2_read.fasta')


# In[ ]:




