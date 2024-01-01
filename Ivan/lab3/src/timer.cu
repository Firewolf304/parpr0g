#include "../include/includes.h"

class timer {
    timespec t1, t2, td;
public:
    __syscall_slong_t NS_PER_SECOND = 1000;
    enum types {
        START = 0,
        FINISH = 1
    };
    string sub_timespec()
    {
        this->td.tv_nsec = this->t2.tv_nsec - this->t1.tv_nsec;
        this->td.tv_sec  = this->t2.tv_sec - this->t1.tv_sec;
        if (this->td.tv_sec > 0 && this->td.tv_nsec < 0)
        {
            this->td.tv_nsec += this->NS_PER_SECOND;
            this->td.tv_sec--;
        }
        else if (this->td.tv_sec < 0 && this->td.tv_nsec > 0)
        {
            this->td.tv_nsec -= this->NS_PER_SECOND;
            this->td.tv_sec++;
        }
        return std::to_string(this->td.tv_sec) + "." + std::to_string(this->td.tv_nsec);
    }
    void operator<< ( types type ) {
        switch(type) {
            case START: {
                clock_gettime(CLOCK_REALTIME, &this->t1);
            }
            case FINISH: {
                clock_gettime(CLOCK_REALTIME, &this->t2);
            }
        }
    }
};