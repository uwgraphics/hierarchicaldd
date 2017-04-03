#ifndef __COMMAND_H__
#define __COMMAND_H__
#include <condition_variable>
#include <exception>
#include <mutex>
#include <queue>
#include <thread>
#include <chrono>
#include <iostream>
namespace Domain_Decomposition{
enum COMMAND{
    IDLE,
    INIT=1,
    SEND_FLAGS,
    SET_CHANNELS,
    STEP1,
    STEP2,
    V_CYCLE,
    SYNC,
    RESTRICTION_UNIT_TEST,
    PROLONGATION_UNIT_TEST,
    BOTTOM_SMOOTH_UNIT_TEST,
    V_CYCLE_TEST_BENCH,
    INTERFACE_COLLECT_DISTRIBUTE_UNIT_TEST
};
enum TAG{
    COMMAND_TAG=1,
    DATA_TAG,
    TERMINATE_TAG
};

template<typename T>
struct Thread_Safe_Queue{
    std::queue<T> queue;
    std::mutex mutex;
    void push(T entry){
        std::unique_lock<std::mutex> lock(mutex);
        queue.push(entry);
    }
    T pop(bool wait=true){
        bool poped = false;
        T entry;
        do{
            std::this_thread::sleep_for(std::chrono::milliseconds(10));
            std::unique_lock<std::mutex> lock(mutex);
            if(queue.empty()) {continue;}
            entry = queue.front();
            queue.pop();
            poped = true;
            break;
        }while(wait);
        if(!poped) throw std::runtime_error("Nothing Poped!");
        return entry;
    }
};
};
#endif
