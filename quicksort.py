
def qsort(a, lo, hi):
    print(a)
    if(lo >= hi):
        return
    p = a[(lo + hi) // 2]       # pivot, any a[] except a[hi]
    print(p)
    i = lo - 1
    j = hi + 1
    while(1):
        while(1):               # while(a[++i] < p)
            i += 1
            if(a[i] >= p):
                break
        while(1):               # while(a[--j] < p)
            j -= 1
            if(a[j] <= p):
                break
        if(i >= j):
            break
        a[i],a[j] = a[j],a[i]
        print(a)
    print(j)
    qsort(a, lo, j)
    qsort(a, j+1, hi)
    
arr = ['D', 'B', 'A', 'C', 'F', 'E', 'G']
# to make best case, have a sorted array use middle as pivot element, do partition, append pivot to array in the order it is seen.
# [a, b, c, d, e, f, g]
# d is pivot. left = [a, b, c], right = [e, f, g] -> array = [d]
# left: b is pivot. left = [a], right = [c] -> array = [d, b] 
# 		-> left([a]) and right([c]) are only 1 element so pivot is themself 
# 		-> array = [d, b, a, c]
# right: f is pivot. left = [e], right = [g] -> array = [d, f] 
# 		-> left([e]) and right([g]) are only 1 element so pivot is themself 
# 		-> array = [d, b, a, c, f, e, g]
qsort(arr, 0, 6)