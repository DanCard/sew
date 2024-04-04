#ifndef THREADPOOL_H
#define THREADPOOL_H

#include <vector>
#include <queue>
#include <thread>
#include <mutex>
#include <condition_variable>
#include <future>
#include <functional>
#include <stdexcept>

class ThreadPool {
public:
  ThreadPool(int);
  template<class F, class... Args>
  auto add_task(F&& f, Args&&... args) 
      -> std::future<typename std::result_of<F(Args...)>::type>;
  ~ThreadPool();

private:
  // Need to keep track of threads so we can join them
  std::vector<std::thread> threads_;
  // The task queue
  std::queue<std::function<void()>> task_q;

  // Synchronization
  std::mutex queue_mutex;
  std::condition_variable condition_var;
  bool stop;
};

// The constructor creates threads
inline ThreadPool::ThreadPool(int threads) : stop(false) {
  for(int i = 0;i<threads;++i)
    threads_.emplace_back([this] {
      for(;;) {
        std::function<void()> task; {
          std::unique_lock<std::mutex> lock(this->queue_mutex);
          this->condition_var.wait(lock,
              [this]{ return this->stop || !this->task_q.empty(); });
          if(this->stop && this->task_q.empty())
              return;
          task = std::move(this->task_q.front());
          this->task_q.pop();
        }
        task();
      }
    }
  );
}

// Add new work item to the pool
template<class F, class... Args>
auto ThreadPool::add_task(F&& f, Args&&... args) 
  -> std::future<typename std::result_of<F(Args...)>::type> {
  using return_type = typename std::result_of<F(Args...)>::type;

  auto task = std::make_shared< std::packaged_task<return_type()> >(
          std::bind(std::forward<F>(f), std::forward<Args>(args)...)
      );
  
  std::future<return_type> future_res = task->get_future(); {
    std::unique_lock<std::mutex> lock(queue_mutex);

    // don't allow adding tasks after stopping the pool
    if(stop)
        throw std::runtime_error("Attempt to add task on stopped ThreadPool");

    task_q.emplace([task](){ (*task)(); });
  }
  condition_var.notify_one();
  return future_res;
}

// The destructor joins all threads
inline ThreadPool::~ThreadPool() {
  {
      std::unique_lock<std::mutex> lock(queue_mutex);
      stop = true;
  }
  condition_var.notify_all();
  for(std::thread &worker: threads_)
      worker.join();
}

#endif

/*
https://www.youtube.com/watch?v=6re5U82KwbY

Thread Pool in C++ is used to manage and efficiently resort to a group (or pool) of threads.
Instead of creating threads again and again for each task and then later destroying them,
what a thread pool does is it maintains a set of pre-created threads.
These threads can be reused again to do many tasks concurrently.
By using this approach we can minimize the overhead that costs us due to the creation and destruction of threads.
This makes our application more efficient.

There is no in-built library in C++ that provides the thread pool,
so we need to create the thread pool manually according to our needs.

What is Thread Pool?
A group of worker threads that are established at the program start and stored in a pool to be used
at a later time are called thread pools.
The Thread Pool effectively maintains and allocates existing threads to do several tasks concurrently,
saving time compared to starting a new thread for each activity.
*/