#pragma once

#include <thread>
#include <functional>
#include <mutex>

#include "Types.h"

//To use this in another class, use CRTP : derive that class (the Owner) from Threads as public Threads<Owner>
//The methods in Threads then become available in the Owner, and std::function pointers can be used in Threads without having to pass the this pointer every time.

template <class Owner> 
class Threads 
{

private:

	Owner *owner;

	//threads are launched as non-blocking by default. For the next thread launch set this to true for blocking. 
	//You must know the threadId, and flag is reset back to default after that threadId is launched.
	std::vector<bool> blocking_call;

	std::vector<long> thread_active, thread_running;

	std::mutex thread_mutex;

protected:

	//reserved unique thread ids for various uses
	enum THREAD_ { THREAD_GENERATENEW = 0, THREAD_LOOP, THREAD_GETINPUT, THREAD_HANDLEMESSAGE, THREAD_HANDLEMESSAGE2, THREAD_HANDLEMESSAGE3, THREAD_HANDLECLIENT, THREAD_NETWORK, THREAD_TIMEDCHECK, THREAD_TIMEDREFRESH };

private:

	//make sure everything is correctly set for the required threadId (or generate a new threadId if asked to - default). Return threadId.
	int configure_threadId(int threadId = THREAD_GENERATENEW);

	int run_on_thread(std::function<void()> runThis, int threadId);

protected:	//all methods usable by other objects are set as protected to force usage of Threads only through CRTP

	Threads(void);
	virtual ~Threads() { Stop_All_Threads(); }

	void Stop_All_Threads() { for (int i = 1; i < (int)thread_active.size(); i++) stop_thread(i); }

	void stop_thread(int threadId);

	int set_blocking_thread(int threadId = THREAD_GENERATENEW);

	bool is_thread_running(int threadId) { if (threadId < (int)thread_active.size()) return (thread_active[threadId] != THREAD_GENERATENEW); else return false; }

	//--------------------INFINITE LOOP THREAD
	//infinite while loop thread : launch it with an Owner method, this thread will call the method without delay in a while loop

	//launch new thread (either using a reserved thread id, or generate a new one) and return thread id
	int infinite_loop_launch(void (Owner::*runThisMethod)(), int threadId = THREAD_GENERATENEW);
	template <typename ... PType> int infinite_loop_launch(void (Owner::*runThisMethod)(PType...), PType... param, int threadId = THREAD_GENERATENEW);

	//--------------------DELAYED CALL THREAD
	//delayed call thread : run the Owner method a single time after a fixed time delay

	int delayed_call_launch(void (Owner::*runThisMethod)(), int delay_ms, int threadId = THREAD_GENERATENEW);
	template <typename ... PType> int delayed_call_launch(void (Owner::*runThisMethod)(PType...), PType... param, int delay_ms, int threadId = THREAD_GENERATENEW);

	//--------------------SINGLE CALL THREAD
	//single call thread : run the Owner method a single time without any delay

	int single_call_launch(void (Owner::*runThisMethod)(), int threadId = THREAD_GENERATENEW);
	template <typename ... PType> int single_call_launch(void (Owner::*runThisMethod)(PType...), PType... param, int threadId = THREAD_GENERATENEW);

	//--------------------TIMED CALL THREAD
	//timed call thread : run the Owner method at fixed time intervals for a given amount of time (or indefinitely if time set to zero)

	int timed_call_launch(void (Owner::*runThisMethod)(), int refresh_ms, int total_ms = 0, int threadId = THREAD_GENERATENEW);
	template <typename ... PType> int timed_call_launch(void (Owner::*runThisMethod)(PType...), PType... param, int refresh_ms, int total_ms = 0, int threadId = THREAD_GENERATENEW);
};

///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

//Methods implemented here, not in a separate .cpp file., due to templating. 
//Alternatively you can implement them in a .cpp file and include that .cpp file in other .cpp files that instantiate this class (stops linker errors as the code for that class is then compiled there).
//Yet another way is to move this code to Threads.cpp, exclude it from the build and #include it here - i.e. effectively achieving same effect as first method.
//Since this file is not large, I prefer the former method as there are less includes to worry about.

template <class Owner> Threads<Owner>::Threads(void) 
{
	//use CRTP to get pointer to derived class (Owner)
	owner = static_cast<Owner*>(this);

	blocking_call.push_back(false);
	thread_active.push_back(THREAD_GENERATENEW);
	thread_running.push_back(THREAD_GENERATENEW);
}

template <class Owner> int Threads<Owner>::configure_threadId(int threadId)
{
	// if this threadId has not been used before, make space for its flags
	if (threadId >= (int)thread_active.size()) {

		thread_active.resize(threadId + 1, 0);
		thread_running.resize(threadId + 1, 0);
		blocking_call.resize(threadId + 1, false);
	}

	else if (threadId == THREAD_GENERATENEW) {

		//generate a new thread id : find the first entry in thread_active which is not being used, or extend the std::vector
		int idx = 1;
		for (; idx < (int)thread_active.size(); idx++)
			if (thread_active[idx] == 0) break;

		//no free entries found : generate new one
		if (idx == (int)thread_active.size()) {

			thread_active.push_back(0);
			thread_running.push_back(0);
			blocking_call.push_back(false);
		}

		//new thread id generated
		threadId = idx;
	}

	return threadId;
}

template <class Owner> void Threads<Owner>::stop_thread(int threadId)
{
	//if not running (or not a valid threadId) then return
	if (!is_thread_running(threadId)) return;

	//stop the thread
	InterlockedDecrement(&thread_running[threadId]);

	//wait for the loop to actually stop - the Sleep is needed to give it a chance to write the value.
	while (thread_active[threadId]) { Sleep(1); }
}

template <class Owner> int Threads<Owner>::set_blocking_thread(int threadId)
{
	thread_mutex.lock();

	threadId = configure_threadId(threadId);

	blocking_call[threadId] = true;

	thread_mutex.unlock();

	return threadId;
}

template <class Owner> int Threads<Owner>::run_on_thread(std::function<void()> runThis, int threadId) 
{
	//std::mutex must be locked before calling this std::function

	//don't start another thread of this type - return false, i.e. couldn't start
	if (thread_active[threadId]) { thread_mutex.unlock(); return THREAD_GENERATENEW; }

	//mark thread_active
	InterlockedIncrement(&thread_active[threadId]);

	//thread will shortly be running
	InterlockedIncrement(&thread_running[threadId]);

	//release the lock now as the call may be blocking, or take a long time to execute - don't want to stop other threads from being created.
	thread_mutex.unlock();

	thread threaded_call(runThis);

	//now launch the infinite loop thread
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
		InterlockedDecrement(&this->thread_active[threadId]);
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
		InterlockedDecrement(&this->thread_active[threadId]);
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

		Sleep(delay_ms);

		CALLFP(this->owner, runThisMethod)();

		//this was set to 1 by the launcher
		InterlockedDecrement(&this->thread_active[threadId]);
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

	std::function<void()> calling_method = [this, param..., delay_ms, threadId, runThisMethod] {

		Sleep(delay_ms);

		CALLFP(this->owner, runThisMethod)(param...);

		//this was set to 1 by the launcher
		InterlockedDecrement(&this->thread_active[threadId]);
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
		InterlockedDecrement(&this->thread_active[threadId]);
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
		InterlockedDecrement(&this->thread_active[threadId]);
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

		DWORD time_elapsed_ms = 0;

		DWORD start_time_ms = GetTickCount();

		while (time_elapsed_ms < (DWORD)total_ms || total_ms == 0) {

			Sleep(refresh_ms);

			CALLFP(this->owner, runThisMethod)();

			time_elapsed_ms = GetTickCount() - start_time_ms;
		}

		//this was set to 1 by the launcher
		InterlockedDecrement(&this->thread_active[threadId]);
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

	std::function<void()> calling_method = [this, param..., refresh_ms, total_ms, threadId, runThisMethod] {

		DWORD time_elapsed_ms = 0;

		DWORD start_time_ms = GetTickCount();

		while (time_elapsed_ms < (DWORD)total_ms || total_ms == 0) {

			Sleep(refresh_ms);

			CALLFP(this->owner, runThisMethod)(param...);

			time_elapsed_ms = GetTickCount() - start_time_ms;
		}

		//this was set to 1 by the launcher
		InterlockedDecrement(&this->thread_active[threadId]);
	};

	return run_on_thread(calling_method, threadId);
}