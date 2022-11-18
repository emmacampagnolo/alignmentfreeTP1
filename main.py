from loading import load_directory
from kmers import stream_kmers, kmer2str, str2kmer
from collections import Counter
import time
import numpy as np
import pandas as pd
import dataframe_image as dfi

def similarity(A, inter, B):
    # --- To complete ---
    return (len(inter)/(len(A)),len(inter)/(len(B)))


def jaccard(A, inter, B):
    # --- To complete ---
    return (len(inter)/(len(A)+len(B)-len(inter)))

def my_method1(file1,file2,k):
    inter=[]
    A=stream_kmers("".join(file1),k)
    A_dict=Counter(A)
    seq="".join(file2)
    B=[]
    kmer=0
    revkmer=0
    mask=(1<<((k-1)*2))-1
    revmask=(1<<(k*2))-1-3
    scores={"A":0,"C":1,"T":2,"G":3}
    rscores={"A":2,"C":3,"T":0,"G":1}
    for i in range(k-1):
        kmer=kmer<<2
        revkmer=revkmer>>2
        kmer+=str2kmer(seq[i],scores)
        revkmer+=str2kmer(seq[i],rscores)<<(k-1)*2
    for nucl in seq[k-1:]:
        kmer=kmer&mask
        #revkmer=revkmer&revmask
        kmer=kmer<<2
        revkmer=revkmer>>2
        kmer+=str2kmer(nucl,scores)
        revkmer+=str2kmer(nucl,rscores)<<(k-1)*2
        Kmer=min(kmer,revkmer)
        B.append(Kmer)
        if Kmer in A_dict:
            inter.append(Kmer)
            if A_dict[Kmer]==1:
                del A_dict[Kmer]
            else:
                A_dict[Kmer]=A_dict[Kmer]-1
    return A, inter, B



def my_method2(file1,file2,k):
    A=sorted(stream_kmers("".join(file1),k))
    B=sorted(stream_kmers("".join(file2),k))
    i=0
    j=0
    inter=[]
    while i<len(A) and j<len(B):
        if A[i]==B[j]:
            inter.append(A[i])
            i+=1
            j+=1
        elif A[i]<B[j]:
            i+=1
        else:
            j+=1
    return A,inter,B
        
    
    


if __name__ == "__main__":
    # Load all the files in a dictionary
    files = load_directory("data")
    k=21
    filenames = list(files.keys())
    adj_matrix=np.zeros((len(filenames),len(filenames)))
    print("Method 1")
    for i in range(len(files)):
        for j in range(i+1, len(files)):
            start_t=time.time()
            A,inter,B=my_method1(files[filenames[i]], files[filenames[j]], k)
            end_t=time.time()
            Jacc=jaccard(A, inter, B)
            print(filenames[i], filenames[j], Jacc, similarity(A, inter, B))
            print("Time spent:", end_t-start_t)
            adj_matrix[i][j]=1-Jacc
            adj_matrix[j][i]=1-Jacc
    
    print("Method 2")
    for i in range(len(files)):
        for j in range(i+1, len(files)):
            start_t=time.time()
            A,inter,B=my_method2(files[filenames[i]], files[filenames[j]], k)
            end_t=time.time()
            print(filenames[i], filenames[j], jaccard(A, inter, B), similarity(A, inter, B))
            print("Time spent:", end_t-start_t)
    
    np.savetxt("adjency_matrix.txt",adj_matrix)
    dfi.export(pd.DataFrame(data=adj_matrix, index=filenames, columns=filenames), "results.png")