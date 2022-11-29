from loading import load_directory
from kmers import stream_kmers, kmer2str, str2kmer
from collections import Counter
import time
import numpy as np
import pandas as pd
import dataframe_image as dfi
import heapq

def similarity(A, inter, B):
    # --- To complete ---
    return (inter/(len(A)),inter/(len(B)))


def jaccard(A, inter, B):
    # --- To complete ---
    return (inter/(len(A)+len(B)-inter))

def my_method1(file1,file2,k):
    inter=[]
    A=stream_kmers("".join(file1),k)
    A_dict=Counter(A)
    seq="".join(file2)
    B=[]
    kmer=0
    revkmer=0
    mask=(1<<((k-1)*2))-1
    scores={"A":0,"C":1,"T":2,"G":3}
    rscores={"A":2,"C":3,"T":0,"G":1}
    for i in range(k-1):
        kmer=kmer<<2
        revkmer=revkmer>>2
        kmer+=str2kmer(seq[i],scores)
        revkmer+=str2kmer(seq[i],rscores)<<(k-1)*2
    for nucl in seq[k-1:]:
        kmer=kmer&mask
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
    A=stream_kmers("".join(file1),k)
    A.sort()
    B=stream_kmers("".join(file2),k)
    B.sort()
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
        

def make_sketch_1(s,seq,k):
    L=[]
    seq="".join(seq)
    for kmer,revkmer in stream_kmers(seq, k):
        Kmer=min(kmer,revkmer)
        Kmer=xorshift64(Kmer)
        if len(L)<s:
            L.append(Kmer)
        else :
            Max=max(L)
            if Kmer<Max:
                imax=L.index(Max)
                L[imax]=Kmer
    L.sort()
    return L

def make_sketch_2(s,seq,k):
    L=[]
    seq="".join(seq)
    for kmer,revkmer in stream_kmers(seq,k):
        Kmer=min(kmer,revkmer)
        Kmer=xorshift64(Kmer)
        if len(L)<s:
            L.append(-Kmer)
            if len(L)==s:
                heapq.heapify(L)
        else:
            Max=L[0]
            if Kmer<-Max:
                heapq.heappushpop(L,-Kmer)
    L=[-x for x in L]
    L.sort()
    return L


def make_sketch_3(s,seq,k):
    L=[float('inf') for i in range(s)]
    seq="".join(seq)
    for kmer,revkmer in stream_kmers(seq,k):
        Kmer=min(kmer,revkmer)
        Kmer=xorshift64(Kmer)
        id_kmer=Kmer%s
        if Kmer<L[id_kmer]:
            L[id_kmer]=Kmer
    L.sort()
    return L
    

def xorshift64(x):
    x ^= x << 13
    x ^= x >> 7
    x ^= x << 17
    return x

def compare_sorted_lists(A,B):
    inter=0
    i=0
    j=0
    while i<len(A) and j<len(B):
        if A[i]==B[j]:
            inter+=1
            i+=1
            j+=1
        elif A[i]<B[j]:
            i+=1
        else:
            j+=1
    return inter

if __name__ == "__main__":
    # Load all the files in a dictionary
    files = load_directory("data")
    k=21
    filenames = list(files.keys())
    s=1000

    for i in range(len(files)):
        for j in range(i+1, len(files)):
            A=make_sketch_1(s, files[filenames[i]], k)
            B=make_sketch_1(s, files[filenames[j]], k)
            inter=compare_sorted_lists(A,B)
            Jacc=jaccard(A, inter, B)
            print(filenames[i], filenames[j], Jacc, similarity(A, inter, B))
            A=make_sketch_2(s, files[filenames[i]], k)
            B=make_sketch_2(s, files[filenames[j]], k)
            inter=compare_sorted_lists(A,B)
            Jacc=jaccard(A, inter, B)
            print(filenames[i], filenames[j], Jacc, similarity(A, inter, B))
            A=make_sketch_3(s, files[filenames[i]], k)
            B=make_sketch_3(s, files[filenames[j]], k)
            inter=compare_sorted_lists(A,B)
            Jacc=jaccard(A, inter, B)
            print(filenames[i], filenames[j], Jacc, similarity(A, inter, B))
            
            
    """
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
    """
 