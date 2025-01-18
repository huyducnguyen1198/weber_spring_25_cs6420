def create_list_of_numbers(n):
	import random
	return [random.randint(-n, n) for _ in range(n)]
	#return [-i for i in range(n//2,0, -1)] + [i for i in range(1, n//2+1)]


def generate_list_unique_sum(size):
    result = []
    current = 1  # Start from the smallest positive integer
    
    while len(result) < size:
        # Check if the current value satisfies the conditions
        valid = True
        sums = {a + b for i, a in enumerate(result) for b in result[i+1:]}  # Compute all pairwise sums
        
        # Ensure no existing sums or potential sums include `current`
        if current in sums:
            valid = False
        else:
            for existing in result:
                if current + existing in result or current + existing in sums:
                    valid = False
                    break
        
        # Add the value if it's valid
        if valid:
            result.append(current)
        
        # Move to the next candidate value
        current += 1
    
    return result

def four_sum_brute(arr):
	total_pairs = 0
	all_pairs = []
	for i in range(len(arr)):
		for j in range(i+1, len(arr)):
			for k in range(j+1, len(arr)):
				for l in range(k+1, len(arr)):
					if arr[i] + arr[j] + arr[k] + arr[l] == 0:
						total_pairs += 1
						all_pairs.append(sorted([arr[i], arr[j], arr[k], arr[l]]))
	
	return total_pairs

def three_sum(arr):
	sum_of_three = {}
	for i in range(len(arr)):
		for j in range(i+1, len(arr)):
			for k in range(j+1, len(arr)):
				suum = arr[i] + arr[j] + arr[k]
				if suum in sum_of_three:
					sum_of_three[suum].append((arr[i], arr[j], arr[k]))
				else:	
					sum_of_three[suum] = [(arr[i], arr[j], arr[k])]
	return sum_of_three
def four_sum_ver2(arr):
	sum_of_three = three_sum(arr)
	sum_of_four = {}
	for i in range(len(arr)):
		if -arr[i] in sum_of_three:
			for triplet in sum_of_three[-arr[i]]:
				if arr[i] not in triplet:
					fourlet = tuple(sorted([arr[i], triplet[0], triplet[1], triplet[2]]))
					if fourlet not in sum_of_four:
						sum_of_four[fourlet] = 1
	
	
	return len(sum_of_four)

def four_sum(arr):
	sum_of_twos = {}
	for i in range(len(arr)):
		for j in range(i+1, len(arr)):
			suum = arr[i] + arr[j]
			if suum in sum_of_twos:
				sum_of_twos[suum].append((arr[i], arr[j]))
			else:    
				sum_of_twos[suum] = [(arr[i], arr[j])]
	
	sum_of_four = {}
	for sum_two in sum_of_twos:
		if -sum_two in sum_of_twos:
			for pair in sum_of_twos[-sum_two]:
				val_1 = pair[0]
				val_2 = pair[1]
				for pair_2 in sum_of_twos[sum_two]:
					
					if val_1 not in pair_2 and val_2 not in pair_2:
						fourlet = tuple(sorted([val_1, val_2, pair_2[0], pair_2[1]]))
						if fourlet not in sum_of_four:
							sum_of_four[fourlet] = 1

	return len(sum_of_four)




def combination(arr, n, k):
	results = []
 
	def backtrack(start, path):
		if len(path) == k:
			comb = [arr[i - 1] for i in path]
			results.append(comb)
			return
		for i in range(start, n + 1):
			path.append(i)
			backtrack(i + 1, path)
			path.pop()
   
	backtrack(1, [])
	return results

def combinations_iterative(arr, n, k):
	# Handle edge cases
	if k > n or k < 0:
		return []
	if k == 0:
		return [[]]

	result = []
	combination = list(range(1, k + 1))  # Initial combination

	while True:
		#result.append(combination[:])  # Add a copy of the current combination

		result.append([arr[i-1] for i in combination])
		# Find the first item from the end that can be incremented
		for i in range(k - 1, -1, -1):
			if combination[i] != i + n - k + 1:
				break
		else:
			return result  # No more combinations possible

		# Increment this item
		combination[i] += 1

		# Reset all items after it
		for j in range(i + 1, k):
			combination[j] = combination[j - 1] + 1


def four_sum_ver3(arr):
	comb_of_three = combinations_iterative(arr, len(arr), 3)
	sum_of_four = {}
	for i in range(len(arr)):
		for comb in comb_of_three:
			sum_three = sum(comb)
			if -arr[i] == sum_three:
				if arr[i] not in comb:
					fourlet = tuple(sorted([arr[i]] + comb))
					if fourlet not in sum_of_four:
						sum_of_four[fourlet] = 1
	return len(sum_of_four)

def four_sum_ver4(arr):
	comb_of_four = combinations_iterative(arr, len(arr), 4)
	
	num_sum_zero = 0
	for comb in comb_of_four:
		if sum(comb) == 0:
			num_sum_zero += 1	
	return num_sum_zero



import time
import matplotlib.pyplot as plt
import random

def test_and_plot_four_sum():
	# Sizes and timings
	sizes = list(range(10,100, 1))  # Interval of 200
	optimized_times = []
	brute_force_times = []
	optimized_times_ver2 = []
	optimized_times_ver3 = []
	optimized_times_ver4 = []
	for size in sizes:
		arr = list(set(generate_list_unique_sum(size)))

		# Measure optimized time
		start_time = time.time()
		num_pairs_optimized = four_sum(arr)
		optimized_times.append(time.time() - start_time)

		# Measure brute force time
		start_time = time.time()
		num_pairs_brute_force = four_sum_brute(arr)
		brute_force_times.append(time.time() - start_time)

		# Measure optimized time ver2
		start_time = time.time()
		num_pairs_optimized_ver2 = four_sum_ver2(arr)
		optimized_times_ver2.append(time.time() - start_time)
		'''
		# Measure optimized time ver3
		start_time = time.time()
		num_pairs_optimized_ver3 = four_sum_ver3(arr)
		optimized_times_ver3.append(time.time() - start_time)

		# Measure optimized time ver4
		start_time = time.time()
		num_pairs_optimized_ver4 = four_sum_ver4(arr)
		optimized_times_ver4.append(time.time() - start_time)'''
	

		print(f"Size: {size}, Brute Force: {num_pairs_brute_force}, Optimized: {num_pairs_optimized},Optimized Ver2: {num_pairs_optimized_ver2}")#, Optimized Ver3: {num_pairs_optimized_ver3}, Optimized Ver4: {num_pairs_optimized_ver4}")
	# Plotting
	plt.plot(sizes, optimized_times, label="Optimized")
	plt.plot(sizes, brute_force_times, label="Brute Force")
	plt.plot(sizes, optimized_times_ver2, label="Optimized Ver2")
	#plt.plot(sizes, optimized_times_ver3, label="Optimized Ver3")
	#plt.plot(sizes, optimized_times_ver4, label="Optimized Ver4")
	plt.xlabel("Size of array")
	plt.ylabel("Time taken (s)")
	plt.legend()
	plt.show()

test_and_plot_four_sum()
'''
size = 100

arr = list(set(create_list_of_numbers(size)))

# Measure optimized time
start_time = time.time()
num_pairs_optimized = four_sum_ver2(arr)
end_time = time.time()
opt_time = end_time - start_time
# Measure brute force time
start_time = time.time()
num_pairs_brute_force = four_sum_brute(arr)
end_time = time.time()
bf_time = end_time - start_time
print(f"Size: {size}, Optimized: {num_pairs_optimized}, Brute Force: {num_pairs_brute_force}, Optimized Time: {opt_time}, Brute Force Time: {bf_time}")'''
