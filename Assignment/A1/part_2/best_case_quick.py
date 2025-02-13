def best_case_quick_rec(data):
	
	# input is sorted array
	def partition(arr):
		'''start, end inclusive'''
		start = 0
		end = len(arr)
		if end - start <= 1: return arr
		mid = end // 2 if end %2 != 0 else end//2 - 1

		pivot = [arr[mid]]
		left = arr[start:mid]
		right = arr[mid + 1: end]
		combined = partition(left) + pivot +  partition(right)
		# swap mid and start
		combined[start], combined[mid] = combined[mid], combined[start]
		print(combined)
		return combined
	result = partition(data)
	def partition_v2(arr, start, end):
		'''start, end inclusive'''

		if end - start <= 1: return 
		mid = (end - start) // 2 if (end - start) %2 != 0 else (end - start)//2 - 1
		mid = start + mid
		pivot = [arr[mid]]
		left = arr[start:mid]
		right = arr[mid + 1: end]
		partition_v2(arr, start, mid)
		partition_v2(arr, mid + 1, end)
		# swap mid and start
		arr[mid], arr[start] = arr[start], arr[mid]		
		return 
	partition_v2(data, 0, len(data))
	#result = partition(data)
	return ''.join(data)
def best_case_quick_iter(arr):
	i = 0 
	j = len(arr) 
	stack = [(i, j)]
	to_swap = []
	while len(stack) > 0:
		start, end = stack.pop(-1)
		mid = start + ((end - start) // 2 if (end - start ) %2 != 0 else (end - start)//2 - 1)

		if (end - start) > 1 :
			stack.append((start, mid))
			stack.append((mid + 1, end))
			to_swap.append((start,mid))
	print(to_swap)
	for start, mid in reversed(to_swap):
		arr[start], arr[mid] = arr[mid], arr[start]
	
	return ''.join(arr)


targer_str = 'KACDBHGFIJEQMLOPNTSRUV'
targer_arr = [letter for letter in targer_str]
input_str  = 'ABCDEFGHIJKLMNOPQRSTUV'
input_arr = [letter for letter in input_str]


best_case_str = best_case_quick_iter(input_arr)

print(best_case_str)
print(best_case_str == targer_str)
