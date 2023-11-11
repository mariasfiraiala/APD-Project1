// Author: APD team, except where source was noted

#include "helpers.h"
#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <string.h>
#include <sys/param.h>
#include <pthread.h>

#define CONTOUR_CONFIG_COUNT    16
#define FILENAME_MAX_SIZE       50
#define STEP                    8
#define SIGMA                   200
#define RESCALE_X               2048
#define RESCALE_Y               2048

#define CLAMP(v, min, max) if(v < min) { v = min; } else if(v > max) { v = max; }

// Creates a map between the binary configuration (e.g. 0110_2) and the corresponding pixels
// that need to be set on the output image. An array is used for this map since the keys are
// binary numbers in 0-15. Contour images are located in the './contours' directory.
ppm_image **init_contour_map() {
    ppm_image **map = (ppm_image **)malloc(CONTOUR_CONFIG_COUNT * sizeof(ppm_image *));
    if (!map) {
        fprintf(stderr, "Unable to allocate memory\n");
        exit(1);
    }

    for (int i = 0; i < CONTOUR_CONFIG_COUNT; i++) {
        char filename[FILENAME_MAX_SIZE];
        sprintf(filename, "./contours/%d.ppm", i);
        map[i] = read_ppm(filename);
    }

    return map;
}

// Updates a particular section of an image with the corresponding contour pixels.
// Used to create the complete contour image.
void update_image(ppm_image *image, ppm_image *contour, int x, int y) {
    for (int i = 0; i < contour->x; i++) {
        for (int j = 0; j < contour->y; j++) {
            int contour_pixel_index = contour->x * i + j;
            int image_pixel_index = (x + i) * image->y + y + j;

            image->data[image_pixel_index].red = contour->data[contour_pixel_index].red;
            image->data[image_pixel_index].green = contour->data[contour_pixel_index].green;
            image->data[image_pixel_index].blue = contour->data[contour_pixel_index].blue;
        }
    }
}

// Corresponds to step 2 of the marching squares algorithm, which focuses on identifying the
// type of contour which corresponds to each subgrid. It determines the binary value of each
// sample fragment of the original image and replaces the pixels in the original image with
// the pixels of the corresponding contour image accordingly.
void march(ppm_image *image, unsigned char **grid, ppm_image **contour_map, int step_x, int step_y) {
    int p = image->x / step_x;
    int q = image->y / step_y;

    for (int i = 0; i < p; i++) {
        for (int j = 0; j < q; j++) {
            unsigned char k = 8 * grid[i][j] + 4 * grid[i][j + 1] + 2 * grid[i + 1][j + 1] + 1 * grid[i + 1][j];
            update_image(image, contour_map[k], i * step_x, j * step_y);
        }
    }
}

// Calls `free` method on the utilized resources.
void free_resources(ppm_image *image, ppm_image **contour_map, unsigned char **grid, int step_x) {
    for (int i = 0; i < CONTOUR_CONFIG_COUNT; i++) {
        free(contour_map[i]->data);
        free(contour_map[i]);
    }
    free(contour_map);

    for (int i = 0; i <= image->x / step_x; i++) {
        free(grid[i]);
    }
    free(grid);

    free(image->data);
    free(image);
}

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

void *thread_func(void *arg) {
    struct thread_argv *argv = (struct thread_argv *)arg;

    // ------------------------------------------- RESCALE ---------------------------------------------------
    // -------------------------------------------------------------------------------------------------------
    // -------------------------------------------------------------------------------------------------------

    // we only rescale downwards
    if (argv->rescale) {
        uint8_t sample[3];

        int start = argv->thread_id * (double)argv->scaled_image->x / argv->num_threads;
        int end = MIN((argv->thread_id + 1) * (double)argv->scaled_image->x / argv->num_threads, argv->scaled_image->x);

        // use bicubic interpolation for scaling
        for (int i = start; i < end; i++) {
            for (int j = 0; j < argv->scaled_image->y; j++) {
                float u = (float)i / (float)(argv->scaled_image->x - 1);
                float v = (float)j / (float)(argv->scaled_image->y - 1);
                sample_bicubic(argv->image, u, v, sample);

                argv->scaled_image->data[i * argv->scaled_image->y + j].red = sample[0];
                argv->scaled_image->data[i * argv->scaled_image->y + j].green = sample[1];
                argv->scaled_image->data[i * argv->scaled_image->y + j].blue = sample[2];
            }
        }
    }

    pthread_barrier_wait(argv->br);
    

    // -------------------------------------------- SAMPLE ---------------------------------------------------
    // -------------------------------------------------------------------------------------------------------
    // -------------------------------------------------------------------------------------------------------

    // // 2. Sample the grid
    // unsigned char **grid = sample_grid(*(argv->scaled_image), argv->step_x, argv->step_y, SIGMA);

    int p = argv->scaled_image->x / argv->step_x;
    int q = argv->scaled_image->y / argv->step_y;

    int start = argv->thread_id * (double)p / argv->num_threads;
    int end = MIN((argv->thread_id + 1) * (double)p / argv->num_threads, p);

    for (int i = start; i < end; i++) {
        for (int j = 0; j < q; j++) {
            ppm_pixel curr_pixel = argv->scaled_image->data[i * argv->step_x * argv->scaled_image->y + j * argv->step_y];

            unsigned char curr_color = (curr_pixel.red + curr_pixel.green + curr_pixel.blue) / 3;

            if (curr_color > argv->sigma) {
                argv->grid[i][j] = 0;
            } else {
                argv->grid[i][j] = 1;
            }
        }
    }

    argv->grid[p][q] = 0;

    pthread_barrier_wait(argv->br);

    // last sample points have no neighbors below / to the right, so we use pixels on the
    // last row / column of the input image for them
    for (int i = start; i < end; i++) {
        ppm_pixel curr_pixel = argv->scaled_image->data[i * argv->step_x * argv->scaled_image->y + argv->scaled_image->x - 1];

        unsigned char curr_color = (curr_pixel.red + curr_pixel.green + curr_pixel.blue) / 3;

        if (curr_color > argv->sigma) {
            argv->grid[i][q] = 0;
        } else {
            argv->grid[i][q] = 1;
        }
    }

    pthread_barrier_wait(argv->br);

    start = argv->thread_id * (double)q / argv->num_threads;
    end = MIN((argv->thread_id + 1) * (double)q / argv->num_threads, q);
    for (int j = start; j < end; j++) {
        ppm_pixel curr_pixel = argv->scaled_image->data[(argv->scaled_image->x - 1) * argv->scaled_image->y + j * argv->step_y];

        unsigned char curr_color = (curr_pixel.red + curr_pixel.green + curr_pixel.blue) / 3;

        if (curr_color > argv->sigma) {
            argv->grid[p][j] = 0;
        } else {
            argv->grid[p][j] = 1;
        }
    }

    pthread_barrier_wait(argv->br);

    // // 3. March the squares
    // march(*(argv->scaled_image), grid, argv->contour_map, argv->step_x, argv->step_y);

    pthread_exit(NULL);
}

int main(int argc, char *argv[]) {
    if (argc < 4) {
        fprintf(stderr, "Usage: ./tema1 <in_file> <out_file> <P>\n");
        return 1;
    }

    ppm_image *image = read_ppm(argv[1]);
    int step_x = STEP;
    int step_y = STEP;

    // 0. Initialize contour map, scaled image and grid
    // ---------------------------------------------- INIT ---------------------------------------------------
    // -------------------------------------------------------------------------------------------------------
    // -------------------------------------------------------------------------------------------------------
    ppm_image **contour_map = init_contour_map();
    ppm_image *scaled_image;
    int rescale = 0;
    if (image->x <= RESCALE_X && image->y <= RESCALE_Y) {
        scaled_image = image;
    } else {
        scaled_image = (ppm_image *)malloc(sizeof(ppm_image));
        if (!scaled_image) {
            fprintf(stderr, "Unable to allocate memory\n");
            exit(1);
        }

        scaled_image->data = (ppm_pixel*)malloc(RESCALE_X * RESCALE_Y * sizeof(ppm_pixel));
        if (!scaled_image->data) {
            fprintf(stderr, "Unable to allocate memory\n");
            exit(1);
        }

        scaled_image->x = RESCALE_X;
        scaled_image->y = RESCALE_Y;
        rescale = 1;
    }

    int p = scaled_image->x / step_x;
    int q = scaled_image->y / step_y;

    unsigned char **grid = (unsigned char **)malloc((p + 1) * sizeof(unsigned char*));
    if (!grid) {
        fprintf(stderr, "Unable to allocate memory\n");
        exit(1);
    }

    for (int i = 0; i <= p; i++) {
        grid[i] = (unsigned char *)malloc((q + 1) * sizeof(unsigned char));
        if (!grid[i]) {
            fprintf(stderr, "Unable to allocate memory\n");
            exit(1);
        }
    }

    // ---------------------------------------------- INIT ---------------------------------------------------
    // -------------------------------------------------------------------------------------------------------
    // -------------------------------------------------------------------------------------------------------

    int num_threads = atoi(argv[3]);
    pthread_t threads[num_threads];
    int rc;
    void *status;
    struct thread_argv arg[num_threads];
    pthread_barrier_t barrier;
    pthread_barrier_init(&barrier, NULL, num_threads);

    for (int i = 0; i < num_threads; ++i) {
        arg[i].num_threads = num_threads;
        arg[i].thread_id = i;

        arg[i].image = image;
        arg[i].contour_map = contour_map;
        arg[i].scaled_image = scaled_image;
        arg[i].step_x = step_x;
        arg[i].step_y = step_y;
        arg[i].br = &barrier;
        arg[i].rescale = rescale;
        arg[i].sigma = SIGMA;
        arg[i].grid = grid;
    }

    for (int i = 0; i < num_threads; ++i) {
        rc = pthread_create(&threads[i], NULL, thread_func, &arg[i]);

        if (rc) {
            printf("Couldn't create thread %d\n", i);
            exit(-1);
        }
    }

    for (int i = 0; i < num_threads; ++i) {
        rc = pthread_join(threads[i], &status);

        if (rc) {
            printf("Couldn't join thread %d\n", i);
            exit(-1);
        }
    }
    pthread_barrier_destroy(&barrier);

    // 3. March the squares
    march(scaled_image, grid, contour_map, step_x, step_y);

    // 4. Write output
    write_ppm(scaled_image, argv[2]);

    free_resources(scaled_image, contour_map, grid, step_x);

    return 0;
}
