Copyright 2023 Maria Sfiraiala (maria.sfiraiala@stud.acs.upb.ro)

# Multithreaded Marching Squares - Project1

## Description

The project aims to parallelize Marching Squares, a graphical algorithm used for getting pictures outlines, in our case, mainly focused on the `PPM` format.
We are also basing our approach on the `Pthreads` library.

### The Thread Function

The implementation heavily relies on syncronizing the threads using barriers due to the quite serial nature of the algorithm.
As a result, the thread function is made up from multiple stages, fenced by barriers:

1. **Rescale** the image, if needed.

   Sometimes, the image is just too big, so, we'll have to scale it down to the maximum size of 2048 x 2048.
   The scaling algorithm means having two imbricated `for`s that calculate the new values of the pixels.
   We parallelized these by the outer indeces values, using the well-known formula for  intervals:

   ```C
   start = thread_id * (double)N / P;
   end = min((thread_id + 1) * (double)N / P, N);
   ```

1. **Sample** the scaled image into a grid.

   We have three stages for this part of the algorithm, due to the fact that we need three `for` blocks: one for the "meat", the center of the image, and two for the last row and column, which don't have enough sample pixels neighbouring them, so we use exactly these lines in constructing the last points.
   Every `for` is guarded by a barrier and is parallelized using the index interval formula mentioned earlier.

1. **March** through image and complete the contour.

   Finally, the last step to complete is getting the binary code of every grid and updating the resulting image with the pieces coresponding to the binary representation of the zone.
   Another imbricated for parallelized via the outer loop and the interval formula.

### The Thread Arguments

Due to the fact that almost all computation is crammed into the thread function, the thread arguments are also numerous:

```C
struct thread_argv {
    int num_threads;
    int thread_id;
    ppm_image *scaled_image;
    ppm_image *image;
    ppm_image **contour_map;
    unsigned char **grid;
    int step_x;
    int step_y;
    pthread_barrier_t *br;
    int rescale;
    int sigma;
};
```

The barrier, the thread id and the number of threads are essential when it comes to syncronization, as they are used to divide: a) the algorithm into parallelizable sections; b) the interval for the loop indeces.

The initial image, the rescale image, the contour map and the grid represent the preallocated data we are using to come up with the contour image.

The rescale value is used to signalyze whether we are required to rescale the original image.

 > **Note**: All the heap allocations and deallocations are handled by the main thread, after all the other threads have been joined.
 This is due to the fact that even though `malloc()` is a thread-safe function, storing and handling its return value (ergo, the pointer to the heap memory block) is tricky between multiple threads and can easily result in memory corruption/leaks.

## Observations Regarding the Project

One of the easiest projects I have worked on so far!
Really fun using all sorts of tools (ASAN, valgrind, gdb) to get info about the memory I was happily and inocently corrupting.

It was very intersting seeing how threads actually got scheduled:
I had some music on during testing and the firefox thread would abruptly stop and immediately start again as the scheduler was planning both the image processing and the browser.
