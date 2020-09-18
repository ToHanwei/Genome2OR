
from multiprocessing import Queue, Pool, Lock, Process
import multiprocessing

#print(multiprocessing.cpu_count())

def fun(x, Q, L):
    x = x**2
    print('value:', x)
    L.acquire()
    Q.put(x)
    L.release()

def main():
    queue = Queue(20)
    lock = Lock()
    inlist = list(range(20))
    for i in inlist[:20:6]:
        jobs = []
        joblist = inlist[i:i+6]
        for i in joblist:
            p = Process(target=fun, args=(i, queue, lock))
            p.start()
            jobs.append(p)

        for p in jobs:
            p.join()
        
    t = []
    while not queue.empty():
        v = queue.get()
        print(v)
        t.append(v)

    print(t)

main()
