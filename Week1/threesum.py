arr = [30, -40, -20, -10, 40, 0, 10, 5]
def create_list_of_numbers(n):
	import random
	return [random.randint(-1000, 10000) for _ in range(n)]


def three_sum(arr):
	sum_of_two = {}
	for i in range(len(arr)):
		for j in range(i+1, len(arr)):
			suum = arr[i] + arr[j]
			if suum in sum_of_two:
				sum_of_two[suum].append((arr[i], arr[j]))
			else:	
				sum_of_two[suum] = [(arr[i], arr[j])]
	
	total_pairs = 0

	sum_of_three = {}
	for i in range(len(arr)):
		if -arr[i] in sum_of_two:
			for pair in sum_of_two[-arr[i]]:
				if arr[i] not in pair:
					triplet = tuple(sorted([arr[i], pair[0], pair[1]]))
					if triplet not in sum_of_three:
						sum_of_three[triplet] = 1
						total_pairs += 1
	return len(sum_of_three)



def three_sum_brute(arr):
	total_pairs = 0
	all_pairs = []
	for i in range(len(arr)):
		for j in range(i+1, len(arr)):
			for k in range(j+1, len(arr)):
				if arr[i] + arr[j] + arr[k] == 0:
					total_pairs += 1
					all_pairs.append(sorted([arr[i], arr[j], arr[k]]))
	 
	return total_pairs


import time
import matplotlib.pyplot as plt
import random

def test_and_plot_three_sum():
	# Sizes and timings
	sizes = list(range(10, 2001, 200))  # Interval of 200
	optimized_times = []
	brute_force_times = []

	for size in sizes:
		print(f"Testing for size {size}")
		arr = list(set(create_list_of_numbers(size)))

		# Measure optimized time
		start_time = time.time()
		three_sum(arr)
		end_time = time.time()
		optimized_times.append(end_time - start_time)

		# Measure brute force time (smaller sizes only to avoid timeout)
		start_time = time.time()
		three_sum_brute(arr)
		end_time = time.time()
		brute_force_times.append(end_time - start_time)

	print(sizes)
	print(optimized_times)
	print(brute_force_times)
	# Plot results
	plt.figure(figsize=(12, 6))
	plt.plot(sizes, optimized_times, label="Optimized Approach", marker='o')

	# Plot brute force times only for available data
 
	
	brute_sizes = [sizes[i] for i in range(len(sizes)) if brute_force_times[i] is not None]
	brute_times = [time for time in brute_force_times if time is not None]
	plt.plot(brute_sizes, brute_times, label="Brute Force Approach", marker='x' )

	plt.xlabel("Input Size (n)")
	plt.ylabel("Time Taken (seconds)")
	plt.title("Three Sum Performance: Optimized vs Brute Force")
	plt.legend()
	plt.grid()
	plt.show()


test_and_plot_three_sum()