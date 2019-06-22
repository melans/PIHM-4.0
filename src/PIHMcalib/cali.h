#ifndef _CALI_H
#define _CALI_H
#include <unistd.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <stdio.h>       /* standard I/O routines.               */
#include <sys/types.h>   /* various type definitions.            */
#include <sys/ipc.h>     /* general SysV IPC structures          */
#include <sys/shm.h>     /* semaphore functions and structs.     */
#include <sys/sem.h>     /* shared memory functions and structs. */
#include <unistd.h>      /* fork(), etc.                         */
#include <sys/wait.h>    /* wait(), etc.                         */
#include <time.h>        /* nanosleep(), etc.                    */
#include <stdlib.h>      /* rand(), etc.                         */

#define SEM_ID    250    /* ID for the semaphore.               */

#include <stdio.h>
#include <errno.h>
#include <fcntl.h>
#include <signal.h>

#define DIM 12
union 
semun {        
    int              val;    /* Value for SETVAL */
    struct semid_ds *buf;    /* Buffer for IPC_STAT, IPC_SET */
    unsigned short  *array;  /* Array for GETALL, SETALL */
    struct seminfo  *__buf;  /* Buffer for IPC_INFO
                                (Linux specific) */
};  

void
random_delay()
{   
    static int initialized = 0;
    int random_num; 
    struct timespec delay;            /* used for wasting time. */
    
    if (!initialized) {
        srand(time(NULL));
        initialized = 1;
    }
    
    random_num = rand() % 10;
    delay.tv_sec = 0;
    delay.tv_nsec = 10*random_num;
    nanosleep(&delay, NULL);
}
void
sem_lock(int sem_set_id)
{   
    /* structure for semaphore operations.   */
    struct sembuf sem_op;
    
    /* wait on the semaphore, unless it's value is non-negative. */
    sem_op.sem_num = 0;
    sem_op.sem_op = -1;
    sem_op.sem_flg = 0;
    semop(sem_set_id, &sem_op, 1);
}
void
sem_unlock(int sem_set_id)
{
    /* structure for semaphore operations.   */
    struct sembuf sem_op;

    /* signal the semaphore - increase its value by one. */
    sem_op.sem_num = 0;
    sem_op.sem_op = 1;   /* <-- Comment 3 */
    sem_op.sem_flg = 0;
    semop(sem_set_id, &sem_op, 1);
}

#endif
