#pragma once

#include <thread>
#include <functional>
#include <mutex>
#include <chrono>
#include <atomic>

#include "Introspection_base.h"

//To use this in another class, use CRTP : derive that class (the Owner) from Threads as public Threads<Owner>
//The methods in Threads then become available in the Owner, and std::function pointers can be used in Threads without having to pass the this pointer every time.

template <class Owner>
class Threads
{

private:

	//Maximum number of concurrent threads possible
	//It's possible to redesign this class so number of maximum threads is dynamically adjusted, but currently not worth the extra complication.
	//The problem is if std::atomic is stored in a vector, since std::atomic is not copy-constructable the std::vector cannot be resized dynamically.
	//It's possible however to have a wrapper for std::atomic which defines a copy constructor and assignment operator, and have this stored in the vector instead
	//You still have the problem of ensuring any resizing is done atomically. 
	//The much simpler solution is to define a maximum number of threads.
	//In the future this could be redesigned, but for now this is fine, the maximum number of threads below is extremely unlikely to be exceeded in any reasonable use case and doesn't take up any significant memory.
	//Also we cannot use a std::vector which is statically sized to THREAD_MAXIMUM, since this is a templated class. 
	//Instead just use good-old pointers and new/delete operators : use RAII, pretty straight-forward to avoid bugs or memory leaks.
	const int THREAD_MAXIMUM = 100;

	Owner *owner;

	//threads are launched as non-blocking by default. For the next thread launch set this to true for blocking. 
	//You must know the threadId, and flag is reset back to default after that threadId is launched.
	bool* blocking_call = nullptr;

	std::atomic<bool>* thread_active = nullptr;
	std::atomic<bool>* thread_running = nullptr;

	std::mutex thread_mutex;

protected:

	//reserved unique thread ids for various uses
	enum THREAD_ { THREAD_ERROR = -1, THREAD_GENERATENEW = 0, THREAD_LOOP, THREAD_LOOP_STOP, THREAD_GETINPUT, THREAD_HANDLEMESSAGE, THREAD_HANDLEMESSAGE2, THREAD_HANDLEMESSAGE3, THREAD_PYTHONSCRIPT, THREAD_DISKACCESS, THREAD_DISPLAY, THREAD_HANDLECLIENT, THREAD_NETWORK, THREAD_TIMEDCHECK, THREAD_TIMEDREFRESH };

private:

	//make sure everything is correctly set for the required threadId (or generate a new threadId if asked to - default). Return threadId.
	int configure_threadId(int threadId = THREAD_GENERATENEW);

	int run_on_thread(std::function<void()> runThis, int threadId);

protected:	//all methods usable by other objects are set as protected to force usage of Threads only through CRTP

	Threads(void);
	virtual ~Threads();

	//Note THREAD_GENERATENEW = 0 is reserved and never used to run any threads
	void Stop_All_Threads() { for (int i = 1; i < THREAD_MAXIMUM; i++) stop_thread(i); }

	void stop_thread(int threadId);

	int set_blocking_thread(int threadId = THREAD_GENERATENEW);
	int set_nonblocking_thread(int threadId = THREAD_GENERATENEW);

	bool is_thread_running(int threadId) const { if (threadId < THREAD_MAXIMUM) return thread_active[threadId]; else return false; }

	//--------------------INFINITE LOOP THREAD
	//infinite while loop thread : launch it with an Owner method, this thread will call the method without delay in a while loop

	//launch new thread (either using a reserved thread id, or generate a new one) and return thread id
	int infinite_loop_launch(void(Owner::*runThisMethod)(), int threadId = THREAD_GENERATENEW);
	//aditionally specify a conditioning function to be called once before starting the while loop
	int infinite_loop_launch(void(Owner::*runThisMethod)(), void(Owner::*conditioningMethod)(), int threadId = THREAD_GENERATENEW);

	template <typename ... PType>
	int infinite_loop_launch(void(Owner::*runThisMethod)(PType...), PType... param, int threadId = THREAD_GENERATENEW);

	template <typename ... PType>
	int infinite_loop_launch(void(Owner::*runThisMethod)(PType...), void(Owner::*conditioningMethod)(), PType... param, int threadId = THREAD_GENERATENEW);

	//--------------------DELAYED CALL THREAD
	//delayed call thread : run the Owner method a single time after a fixed time delay

	int delayed_call_launch(void(Owner::*runThisMethod)(), int delay_ms, int threadId = THREAD_GENERATENEW);

	template <typename ... PType>
	int delayed_call_launch(void(Owner::*runThisMethod)(PType...), PType... param, int delay_ms, int threadId = THREAD_GENERATENEW);

	//--------------------SINGLE CALL THREAD
	//single call thread : run the Owner method a single time without any delay

	int single_call_launch(void(Owner::*runThisMethod)(), int threadId = THREAD_GENERATENEW);

	template <typename ... PType>
	int single_call_launch(void(Owner::*runThisMethod)(PType...), PType... param, int threadId = THREAD_GENERATENEW);

	//--------------------TIMED CALL THREAD
	//timed call thread : run the Owner method at fixed time intervals for a given amount of time (or indefinitely if time set to zero)

	int timed_call_launch(void(Owner::*runThisMethod)(), int refresh_ms, int total_ms = 0, int threadId = THREAD_GENERATENEW);

	template <typename ... PType>
	int timed_call_launch(void(Owner::*runThisMethod)(PType...), PType... param, int refresh_ms, int total_ms = 0, int threadId = THREAD_GENERATENEW);
};

///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

//Methods implemented here, not in a separate .cpp file., due to templating. 
//Alternatively you can implement them in a .cpp file and include that .cpp file in other .cpp files that instantiate this class (stops linker errors as the code for that class is then compiled there).
//Yet another way is to move this code to Threads.cpp, exclude it from the build and #include it here - i.e. effectively achieving same effect as first method.
//Since this file is not large, I prefer the former method as there are less includes to worry about.

template <class Owner>
Threads<Owner>::Threads(void)
{
	//use CRTP to get pointer to derived class (Owner)
	owner = static_cast<Owner*>(this);

	blocking_call = new bool[THREAD_MAXIMUM];
	thread_active = new std::atomic<bool>[THREAD_MAXIMUM];
	thread_running = new std::atomic<bool>[THREAD_MAXIMUM];

	for (int idx = 0; idx < THREAD_MAXIMUM; idx++) {

		blocking_call[idx] = false;
		thread_active[idx] = false;
		thread_running[idx] = false;
	}
}

template <class Owner>
Threads<Owner>::~Threads()
{
	Stop_All_Threads();

	delete[] blocking_call;
	delete[] thread_active;
	delete[] thread_running;
}

template <class Owner>
int Threads<Owner>::configure_threadId(int threadId)
{
	if (threadId >= THREAD_MAXIMUM) return THREAD_ERROR;

	else if (threadId == THREAD_GENERATENEW) {

		//find the first non-active thread
		int idx = 1;
		for (; idx < THREAD_MAXIMUM; idx++) {

			if (!thread_active[idx]) break;
		}

		//no free entries found : problem! Will just have to wait until something becomes available, but highly unlikely this will ever happen in practice.
		if (idx == THREAD_MAXIMUM) return THREAD_ERROR;

		//new thread id generated
		threadId = idx;
	}

	return threadId;
}

template <class Owner>
void Threads<Owner>::stop_thread(int threadId)
{
	//if not running (or not a valid threadId) then return
	if (!is_thread_running(threadId)) return;

	while (thread_active[threadId]) {

		//stop the thread
		std::atomic_store(&thread_running[threadId], false);

		//wait for the loop to actually stop - the sleep is needed to give it a chance to write the value.
		std::this_thread::sleep_for(std::chrono::milliseconds(1));
	}
}

template <class Owner>
int Threads<Owner>::set_blocking_thread(int threadId)
{
	thread_mutex.lock();

	threadId = configure_threadId(threadId);

	if (threadId != THREAD_ERROR) blocking_call[threadId] = true;

	thread_mutex.unlock();

	return threadId;
}

template <class Owner>
int Threads<Owner>::set_nonblocking_thread(int threadId)
{
	thread_mutex.lock();

	threadId = configure_threadId(threadId);

	if (threadId != THREAD_ERROR) blocking_call[threadId] = false;

	thread_mutex.unlock();

	return threadId;
}

template <class Owner>
int Threads<Owner>::run_on_thread(std::function<void()> runThis, int threadId)
{
	//std::mutex must be locked before calling this std::function

	if (threadId == THREAD_ERROR || std::atomic_exchange(&thread_active[threadId], true)) { thread_mutex.unlock(); return THREAD_GENERATENEW; }

	//thread will shortly be running
	thread_running[threadId] = true;

	//release the lock now as the call may be blocking, or take a long time to execute - don't want to stop other threads from being created.
	thread_mutex.unlock();

	std::thread threaded_call(runThis);

	//now launch the thread
	if (blocking_call[threadId]) threaded_call.join();
	else threaded_call.detach();

	//reset blocking call flag
	blocking_call[threadId] = false;

	return threadId;
}

///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//INFINITE LOOP

template <class Owner>
int Threads<Owner>::infinite_loop_launch(void (Owner::*runThisMethod)(), int threadId)
{
	//ensure thread-safe access
	thread_mutex.lock();

	threadId = configure_threadId(threadId);

	std::function<void()> calling_method = [this, threadId, runThisMethod] {

		//infinite_loop_running was incremented by the launcher
		while (this->thread_running[threadId]) {

			CALLFP(this->owner, runThisMethod)();
		}

		//infinite_loop_thread_active was incremented by the launcher
		std::atomic_store(&this->thread_active[threadId], false);
	};

	//the std::mutex will be unlocked before this std::function returns
	return run_on_thread(calling_method, threadId);
}

template <class Owner>
int Threads<Owner>::infinite_loop_launch(void(Owner::*runThisMethod)(), void(Owner::*conditioningMethod)(), int threadId)
{
	//ensure thread-safe access
	thread_mutex.lock();

	threadId = configure_threadId(threadId);

	std::function<void()> calling_method = [this, threadId, runThisMethod, conditioningMethod] {

		CALLFP(this->owner, conditioningMethod)();

		//infinite_loop_running was incremented by the launcher
		while (this->thread_running[threadId]) {

			CALLFP(this->owner, runThisMethod)();
		}

		//infinite_loop_thread_active was incremented by the launcher
		std::atomic_store(&this->thread_active[threadId], false);
	};

	//the std::mutex will be unlocked before this std::function returns
	return run_on_thread(calling_method, threadId);
}

template <class Owner>
template <typename ... PType>
int Threads<Owner>::infinite_loop_launch(void (Owner::*runThisMethod)(PType...), PType... param, int threadId)
{
	//ensure thread-safe access
	thread_mutex.lock();

	threadId = configure_threadId(threadId);

	std::function<void()> calling_method = [this, param..., threadId, runThisMethod] {

		//infinite_loop_running was incremented by the launcher
		while (this->thread_running[threadId]) {

			CALLFP(this->owner, runThisMethod)(param...);
		}

		//infinite_loop_thread_active was incremented by the launcher
		std::atomic_store(&this->thread_active[threadId], false);
	};

	return run_on_thread(calling_method, threadId);
}

template <class Owner>
template <typename ... PType>
int Threads<Owner>::infinite_loop_launch(void (Owner::*runThisMethod)(PType...), void(Owner::*conditioningMethod)(), PType... param, int threadId)
{
	//ensure thread-safe access
	thread_mutex.lock();

	threadId = configure_threadId(threadId);

	std::function<void()> calling_method = [this, param..., threadId, runThisMethod, conditioningMethod] {

		CALLFP(this->owner, conditioningMethod)();
	
		//infinite_loop_running was incremented by the launcher
		while (this->thread_running[threadId]) {

			CALLFP(this->owner, runThisMethod)(param...);
		}

		//infinite_loop_thread_active was incremented by the launcher
		std::atomic_store(&this->thread_active[threadId], false);
	};

	return run_on_thread(calling_method, threadId);
}

///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//DELAYED CALL

template <class Owner>
int Threads<Owner>::delayed_call_launch(void (Owner::*runThisMethod)(), int delay_ms, int threadId)
{
	//ensure thread-safe access
	thread_mutex.lock();

	threadId = configure_threadId(threadId);

	std::function<void()> calling_method = [this, delay_ms, threadId, runThisMethod] {

		std::this_thread::sleep_for(std::chrono::milliseconds(delay_ms));

		CALLFP(this->owner, runThisMethod)();

		//this was set to 1 by the launcher
		std::atomic_store(&this->thread_active[threadId], false);
	};

	return run_on_thread(calling_method, threadId);
}

template <class Owner>
template <typename ... PType>
int Threads<Owner>::delayed_call_launch(void (Owner::*runThisMethod)(PType...), PType... param, int delay_ms, int threadId)
{
	//ensure thread-safe access
	thread_mutex.lock();

	threadId = configure_threadId(threadId);

	std::function<void()> calling_method = [this, param..., delay_ms, threadId, runThisMethod]{

		std::this_thread::sleep_for(std::chrono::milliseconds(delay_ms));

		CALLFP(this->owner, runThisMethod)(param...);

		//this was set to 1 by the launcher
		std::atomic_store(&this->thread_active[threadId], false);
	};

	return run_on_thread(calling_method, threadId);
}

///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//SINGLE CALL

//example call, run HandleCommand(command) method only once, in its own asynchronous thread: threaded_call.single_call_launch<std::string>(&Simulation::HandleCommand, command)

template <class Owner>
int Threads<Owner>::single_call_launch(void (Owner::*runThisMethod)(), int threadId)
{
	//ensure thread-safe access
	thread_mutex.lock();

	threadId = configure_threadId(threadId);

	std::function<void()> calling_method = [this, threadId, runThisMethod] {

		CALLFP(this->owner, runThisMethod)();

		//this was set to 1 by the launcher
		std::atomic_store(&this->thread_active[threadId], false);
	};

	return run_on_thread(calling_method, threadId);
}

template <class Owner>
template <typename ... PType>
int Threads<Owner>::single_call_launch(void (Owner::*runThisMethod)(PType...), PType... param, int threadId)
{
	//ensure thread-safe access
	thread_mutex.lock();

	threadId = configure_threadId(threadId);

	std::function<void()> calling_method = [this, param..., threadId, runThisMethod] {

		CALLFP(this->owner, runThisMethod)(param...);

		//this was set to 1 by the launcher
		std::atomic_store(&this->thread_active[threadId], false);
	};

	return run_on_thread(calling_method, threadId);
}

///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//TIMED CALL

template <class Owner>
int Threads<Owner>::timed_call_launch(void (Owner::*runThisMethod)(), int refresh_ms, int total_ms, int threadId)
{
	//ensure thread-safe access
	thread_mutex.lock();

	threadId = configure_threadId(threadId);

	std::function<void()> calling_method = [this, refresh_ms, total_ms, threadId, runThisMethod] {

		int time_elapsed_ms = 0;

		while (time_elapsed_ms < total_ms || total_ms == 0) {

			std::this_thread::sleep_for(std::chrono::milliseconds(refresh_ms));

			CALLFP(this->owner, runThisMethod)();

			time_elapsed_ms += refresh_ms;
		}

		//this was set to 1 by the launcher
		std::atomic_store(&this->thread_active[threadId], false);
	};

	return run_on_thread(calling_method, threadId);
}

template <class Owner>
template <typename ... PType>
int Threads<Owner>::timed_call_launch(void (Owner::*runThisMethod)(PType...), PType... param, int refresh_ms, int total_ms, int threadId)
{
	//ensure thread-safe access
	thread_mutex.lock();

	threadId = configure_threadId(threadId);

	std::function<void()> calling_method = [this, param..., refresh_ms, total_ms, threadId, runThisMethod]{

		int time_elapsed_ms = 0;

	while (time_elapsed_ms < total_ms || total_ms == 0) {

			std::this_thread::sleep_for(std::chrono::milliseconds(refresh_ms));

			CALLFP(this->owner, runThisMethod)(param...);

			time_elapsed_ms += refresh_ms;
		}

	//this was set to 1 by the launcher
	std::atomic_store(&this->thread_active[threadId], false);
	};

	return run_on_thread(calling_method, threadId);
}