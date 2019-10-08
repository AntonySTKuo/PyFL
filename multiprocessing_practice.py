import time
import numpy as np
import concurrent.futures
from tqdm import tqdm

start = time.perf_counter()

# def do_something(seconds):
# 	#print(f'Sleeping {seconds} second...')
# 	time.sleep(seconds)
# 	return 'Done Sleeping !!'

def do_something2(index):
	#print(f'index = {index}')
	time.sleep(0.01)
	return index+1

a = np.zeros(1000)
with concurrent.futures.ProcessPoolExecutor() as executor:
	indexes = range(1000)
	for index, value in tqdm(zip(indexes, executor.map(do_something2, indexes)), total=len(indexes)):
		a[index] = value
print(a)

# with concurrent.futures.ProcessPoolExecutor() as executor:
# 	secs = [10, 9, 8, 7, 6, 5, 4, 3, 2, 1]
# 	results = executor.map(do_something, secs)

# with concurrent.futures.ProcessPoolExecutor() as executor:
# 	secs = [10, 9, 8, 7, 6, 5, 4, 3, 2, 1]
# 	results = [executor.submit(do_something, sec) for sec in secs]
# 	for f in tqdm(concurrent.futures.as_completed(results)):
# 		test = 1

# with concurrent.futures.ProcessPoolExecutor() as executor:
# 	f1 = executor.submit(do_something, 1)
# 	print(f1.result())

# do_something()
# do_something()

finish = time.perf_counter()

print(f'Finished in {round(finish-start, 2)} second(s)')

