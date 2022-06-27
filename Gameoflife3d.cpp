#include <algorithm>
#include <assert.h>
#include <math.h>
#include <tuple>
#include <mpi.h>
#include <omp.h>
#include <stdarg.h>
#include <stdio.h>
#include <string.h>
#include <unistd.h>

#define matrix_owner(index, process, n) (((process) * ((index) + 1) - 1) / (n))

//  node used for representing the matrix
struct node {
    short z; 
    short neighbours; 
    bool dead; 
    bool operator<(const struct node& h) const
    {
        return z < h.z;
    }
};

typedef struct array_struct {
    size_t initial_size;
    void* data;
    size_t used;
    size_t size;
    size_t step;
} dynamic_array;


inline int max_two_pow(int n)
{
    n--;
    n |= n >> 1; n |= n >> 2; n |= n >> 4; n |= n >> 8; n |= n >> 16;
    n++;
    return n;
}

void array_init(dynamic_array* array, size_t initial_size)
{
    array->initial_size = max_two_pow(initial_size);
    array->step = sizeof(struct node);
    array->used = 0;

    array->size = array->initial_size;
    array->data = realloc(NULL, array->initial_size * array->step);
    assert(array->data);
}

void array_resize(dynamic_array* array, size_t new_size)
{
    assert(array);

    array->size = max_two_pow(new_size);
    array->data = realloc(array->data, array->size * array->step);
    assert(array->data);
    
}

inline void array_contract(dynamic_array* array)
{
    if (array->size > array->initial_size && array->used <= array->size / 4) {
        array_resize(array, array->size / 2);
    }
}

void array_insert(dynamic_array* array, struct node* to_insert)
{
    assert(array);

    if (array->used == array->size) {
        array_resize(array, array->size * 2);
    }

    struct node* dest = ((struct node*)array->data) + array->used;
    memcpy(dest, to_insert, array->step);
    array->used++;
}

void array_delete_pos(dynamic_array* array, size_t i)
{
    assert(array);
  
    if (i >= array->used) { 
        printf("Invalid index array: smaller or larger\n");
        return;
    }

    array->used--;

    struct node* dest = ((struct node*)array->data) + i;
    struct node* src = ((struct node*)array->data) + array->used;
    memcpy(dest, src, array->step);

    array_contract(array);
}

void array_free(dynamic_array* array)
{
    assert(array);

    if (array->data) {
        free(array->data);
        array->data = NULL;
        array->used = 0;
        array->size = 0;
    }
}

void print_node(struct node* n)
{
    printf("{z: %hd, num_nei: %hd, dead: %s}\n", n->z, n->neighbours, n->dead ? "true" : "false");
}

void array_print(dynamic_array* array)
{
    assert(array);

    printf("size: %lu\n", array->size);
    printf("used: %lu\n", array->used);
    printf("data:\n");

    size_t i;
    for (i = 0; i < array->used; i++) {
        struct node* ptr = ((struct node*)array->data) + i;
        printf("  ");
        print_node(ptr);
    }

}

int array_find_z(dynamic_array* array, short test_z)
{
    assert(array);

    fflush(stdout);

    for (size_t i = 0; i < array->used; i++) {
        struct node* ptr = ((struct node*)array->data) + i;

        if (ptr->z == test_z) {
            return i;
        }
    }
    return -1;
}

struct matrix_struct {
    short side;
    dynamic_array** data;
};
typedef matrix_struct Matrix;
void matrix_print(Matrix* m);

// initializing matrix
inline Matrix matrix_create(short side)
{
    Matrix m;
    m.side = side;
    m.data = (dynamic_array**)calloc(side * side, sizeof(dynamic_array*));
    return m;
}

// Returns the first row and column of the matrix
inline dynamic_array* get_matrix(Matrix* m, short x, short y)
{
    return m->data[x + (y * m->side)];
}

// Returns the element from the matrix
inline struct node* get_matrix_element(Matrix* m, short x, short y, short z)
{
    dynamic_array* array = get_matrix(m, x, y);
    if (!array) {
        return NULL;
    }

    int pos = array_find_z(array, z);
    if (pos == -1) {
        return NULL;
    } else {
        return ((struct node*)array->data) + pos;
    }
}


// Remove from the matrix
inline void matrix_remove(Matrix* m, short x, short y, short z)
{
    dynamic_array* array = get_matrix(m, x, y);
    if (!array) {
        return;
    }
    int pos = array_find_z(array, z);
    if (pos != -1) {
        array_delete_pos(array, pos);
    }
}

void matrix_free(Matrix* m)
{
    for (int i = 0; i < m->side; i++) {
        for (int j = 0; j < m->side; j++) {
            dynamic_array* array = get_matrix(m, i, j);
            if (!array) {
                continue;
            }
            array_free(array);
            free(array);
            array = NULL;
        }
    }
    free(m->data);
}

// alive cells in the matrix printed
void matrix_print_live(Matrix* m, int from = 0, int to = -1)
{
    short SIZE = m->side;
    if (to == -1) {
        to = SIZE - 1;
    }
    for (short x = from; x <= to; x++) {
        for (short y = 0; y < SIZE; y++) {
            dynamic_array* array = get_matrix(m, x, y);
            if (array != NULL) {
                struct node* ptr = ((struct node*)array->data);
                std::sort(ptr, (ptr + array->used));
                for (size_t k = 0; k < array->used; k++) {
                    printf("%hd %hd %hd\n", x, y, ptr->z);
                    ptr++;
                }
            }
        }
    }
}

#define matrix_init 17
#define length_z 16
#define rows_swap 18

void get_z_lengths(Matrix* m, int x, int* z_lengths)
{

    #pragma omp parallel for shared(m,x,z_lengths)
    // length of each z list
    for (int y = 0; y < m->side; y++) { 
        dynamic_array* array = get_matrix(m, x, y);
        z_lengths[y] = (array == NULL) ? 0 : array->used;
    }
}

void matrix_print(Matrix* m)
{
    short SIZE = m->side;
    for (short i = 0; i < SIZE; i++) {
        for (short j = 0; j < SIZE; j++) {
            dynamic_array* array = get_matrix(m, i, j);
            if (array != NULL) {
                struct node* ptr = ((struct node*)array->data);
                printf("(%hd, %hd): [", i, j);
                for (size_t k = 0; k < array->used; k++) {
                    printf("%hd%c", ptr->z, k == array->used - 1 ? '\x07' : ',');
                    ptr++;
                }
                printf("] (size: %lu; used: %lu)\n", array->size, array->used);
            } else {
                printf("(%hd, %hd): []\n", i, j);
            }
        }
    }
}

inline short val_mod(short val, short mod)
{
    if (val >= mod)
        return val - mod;
    else if (val < 0)
        return val + mod;
    else
        return val;
}

// sending the source , row to send and dest
void init_send_row(Matrix* m, int x, int to)
{
    int SIZE = m->side;
    int z_lengths[SIZE];
    get_z_lengths(m, x, z_lengths);

    MPI_Send(z_lengths, SIZE, MPI_INT, to, length_z, MPI_COMM_WORLD);

   
    for (int y = 0; y < SIZE; y++) {
        dynamic_array* array = get_matrix(m, x, y);
        if (array == NULL || array->used == 0) { 
            continue;
        }
        MPI_Send(array->data, (z_lengths[y] * array->step), MPI_BYTE, to, matrix_init, MPI_COMM_WORLD);
    }
}
#define max(value1, value2) (((value1) > (value2)) ? (value1) : (value2))

dynamic_array* array_make_ptr(size_t initial_size)
{
    dynamic_array* array = (dynamic_array*)malloc(sizeof(dynamic_array));
    array_init(array, max(4, initial_size));
    return array;
}

// Insert a new element in the matrix 
inline void matrix_insert(Matrix* m, short x, short y, short z, bool dead, short num_nei)
{
    struct node new_el = { z, num_nei, dead };

    dynamic_array* array = get_matrix(m, x, y);
    if (array == NULL) {
        array = array_make_ptr(4);
        m->data[x + (y * m->side)] = array;
    }
    array_insert(array, &new_el);
}

void init_recv_row(Matrix* m, int x, int owner = 0)
{
    int SIZE = m->side;
    int z_lengths[SIZE];
    MPI_Recv(z_lengths, SIZE, MPI_INT, owner, length_z, MPI_COMM_WORLD, MPI_STATUS_IGNORE);

    for (int y = 0; y < SIZE; y++) {
        dynamic_array* array = get_matrix(m, x, y);
        if (z_lengths[y] == 0) {
            if (array) {
                array->used = 0;
            }
            continue;
        }

        if (array == NULL) {
            array = array_make_ptr(z_lengths[y]);
            m->data[x + (y * SIZE)] = array;
        } else {
            array_resize(array, z_lengths[y]);
        }

        array->used = z_lengths[y];
        MPI_Recv(array->data, (array->used * array->step), MPI_BYTE, owner, matrix_init, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
    }
}

#define matrix_down(id, process, n) ((id) * (n) / (process))
#define matrix_up(id, process, n) (matrix_down((id) + 1, process, n) - 1)
#define matrix_size(id, process, n) (matrix_up(id, process, n) - matrix_down(id, process, n) + 1)

void swap_rows(Matrix* m, int x_have, int x_want, int to, int tmp_id)
{
    MPI_Request request;
    MPI_Status status;
    int SIZE = m->side;

    // receiving all the lengths 
    int their_z_lengths[SIZE];
    MPI_Irecv(their_z_lengths, SIZE, MPI_INT, to, length_z, MPI_COMM_WORLD, &request);

    // Sending all the lengths
    int my_z_lengths[SIZE];
    get_z_lengths(m, x_have, my_z_lengths);
    MPI_Send(my_z_lengths, SIZE, MPI_INT, to, length_z, MPI_COMM_WORLD);

    // Wait for completion of the recv
    MPI_Wait(&request, &status);

    MPI_Request requests[SIZE];
    MPI_Status statuses[SIZE];
    int request_i = 0; 
   

    // Receive all z lists 
    for (int y = 0; y < SIZE; y++) {
        dynamic_array* array = get_matrix(m, x_want, y);
        if (their_z_lengths[y] == 0) {
            if(array) {
                array->used = 0;
            }
            continue;
        }

        if (array == NULL) {
            array = array_make_ptr(their_z_lengths[y]);
            m->data[x_want + (y * SIZE)] = array;
        } else {
            array_resize(array, their_z_lengths[y]);
        }

        array->used = their_z_lengths[y];
        MPI_Irecv(array->data, (array->used * array->step), MPI_BYTE, to, rows_swap, MPI_COMM_WORLD, &requests[request_i]);
        request_i++;
    }

    // Send each z list in the row fashion
    for (int y = 0; y < SIZE; y++) {
        dynamic_array* array = get_matrix(m, x_have, y);
        if (array == NULL || array->used == 0) { 
            continue;
        }

        MPI_Send(array->data, my_z_lengths[y] * array->step, MPI_BYTE, to, rows_swap, MPI_COMM_WORLD);
    }

    if (request_i != 0) {
        MPI_Waitall(request_i, requests, statuses);
    }
}

int main(int argc, char* argv[])
{
    setvbuf(stdout, NULL, _IONBF, 0);

    int id, process;
    double time_passed, divide_time, run_time, collect_time;

    MPI_Init(&argc, &argv);
    MPI_Barrier(MPI_COMM_WORLD);
    time_passed = -MPI_Wtime();

    MPI_Comm_rank(MPI_COMM_WORLD, &id);
    MPI_Comm_size(MPI_COMM_WORLD, &process);
    omp_set_num_threads(2);

    Matrix m;
    int generations;
    short SIZE;
    char* input_file;
    if (id == 0) {
        if (argc != 3) {
            printf(" Incorrect command line argument. Use below!\n");
            printf("[Usage] %s [input-file] [iterations]\n", argv[0]);
            return -1;
        }
        input_file = argv[1];
        generations = atoi(argv[2]);

        if (generations <= 0) {
            printf(" generations must be more than zero : '%d'\n", generations);
            return -1;
        }

        FILE* fp = fopen(input_file, "r");
        if (fp == NULL) {
            printf("Unable to read the input file.\n");
            perror("[ERROR]");
            return -1;
        }

        if (fscanf(fp, "%hd", &SIZE) == EOF) {
            printf(" Unable to read the size.\n");
            return -1;
        }

        Matrix aux = matrix_create(SIZE);
        short x, y, z;
        while (fscanf(fp, "%hd %hd %hd", &x, &y, &z) != EOF) {
            matrix_insert(&aux, x, y, z, false, -1);
        }
        fclose(fp);

        
        // Send gen and size
        int buf[2];
        buf[0] = (int)SIZE;
        buf[1] = generations;
        MPI_Bcast(buf, 2, MPI_INT, 0, MPI_COMM_WORLD);

        m = matrix_create(SIZE);

        // Send rows
        for (int x = 0; x < SIZE; x++) {
            int owner = matrix_owner(x, process, SIZE);

            if (owner != 0) {
                init_send_row(&aux, x, owner);
            } else {
                for (int y = 0; y < SIZE; y++) {
                    dynamic_array* array = get_matrix(&aux, x, y);
                    if (array == NULL) {
                        continue;
                    }

                    dynamic_array* new_array = array_make_ptr(array->used);
                    m.data[x + (y * m.side)] = new_array;

                    new_array->used = array->used;
                    memcpy(new_array->data, array->data, array->used * new_array->step);
                }
            }
        }

        matrix_free(&aux);

    } else {
        // Receive size and gen
        int buf[2];
        MPI_Bcast(buf, 2, MPI_INT, 0, MPI_COMM_WORLD);
        SIZE = (short)buf[0];
        generations = buf[1];
        // receive rows
        m = matrix_create(SIZE);
        for (int x = matrix_down(id, process, SIZE); x <= matrix_up(id, process, SIZE); x++) {
            init_recv_row(&m, x);
        }
    }

    int high = matrix_up(id, process, SIZE);
    int low = matrix_down(id, process, SIZE);

    // front rows
    int higher = val_mod(high + 1, SIZE);
    int lower = val_mod(low - 1, SIZE);

    int higher_owner = matrix_owner(higher, process, SIZE);
    int lower_owner = matrix_owner(lower, process, SIZE);

    MPI_Barrier(MPI_COMM_WORLD);
    divide_time = time_passed + MPI_Wtime();

    for (int gen = 0; gen < generations; gen++) {
        dynamic_array* array;
        struct node *ptr, *test;
        int32_t j;
        short z, c, b, a, y, x;
        if (id % 2 == 0) {
            swap_rows(&m, high, higher, higher_owner, id);
            swap_rows(&m, low, lower, lower_owner, id);
        } else {
            swap_rows(&m, low, lower, lower_owner, id);
            swap_rows(&m, high, higher, higher_owner, id);
        }
        //game logic
        for (int a1 = 0, i = lower; a1 < matrix_size(id, process, SIZE) + 2; a1++) {
            for (j = 0; j < SIZE; j++) {
                array = get_matrix(&m, i, j);
                if (!array) {
                    continue;
                }

                size_t limit = array->used;
                for (size_t k = 0; k < limit; k++) {
                    ptr = ((struct node*)array->data) + k;

                    if (ptr->dead) {
                        continue;
                    }

                    x = i;
                    y = j;
                    z = ptr->z;
                    ptr->neighbours = 0;

                    c = val_mod(z + 1, SIZE);
                    test = get_matrix_element(&m, x, y, c);
                    if (test) {
                        if (test->dead == true) {
                            test->neighbours++;
                        } else {
                            ptr->neighbours++;
                        }
                    } else {
                        matrix_insert(&m, x, y, c, true, 1);
                        ptr = ((struct node*)array->data) + k;
                    }

                    c = val_mod(z - 1, SIZE);
                    test = get_matrix_element(&m, x, y, c);
                    if (test) {
                        if (test->dead == true) {
                            test->neighbours++;
                        } else {
                            ptr->neighbours++;
                        }
                    } else {
                        matrix_insert(&m, x, y, c, true, 1);
                        ptr = ((struct node*)array->data) + k;
                    }

                    a = val_mod(x + 1, SIZE);
                    test = get_matrix_element(&m, a, y, z);
                    if (test) {
                        if (test->dead == true) {
                            test->neighbours++;
                        } else {
                            ptr->neighbours++;
                        }
                    } else {
                        matrix_insert(&m, a, y, z, true, 1);
                        ptr = ((struct node*)array->data) + k;
                    }

                    a = val_mod(x - 1, SIZE);
                    test = get_matrix_element(&m, a, y, z);
                    if (test) {
                        if (test->dead == true) {
                            test->neighbours++;
                        } else {
                            ptr->neighbours++;
                        }
                    } else {
                        matrix_insert(&m, a, y, z, true, 1);
                        ptr = ((struct node*)array->data) + k;
                    }

                    b = val_mod(y + 1, SIZE);
                    test = get_matrix_element(&m, x, b, z);
                    if (test) {
                        if (test->dead == true) {
                            test->neighbours++;
                        } else {
                            ptr->neighbours++;
                        }
                    } else {
                        matrix_insert(&m, x, b, z, true, 1);
                        ptr = ((struct node*)array->data) + k;
                    }

                    b = val_mod(y - 1, SIZE);
                    test = get_matrix_element(&m, x, b, z);
                    if (test) {
                        if (test->dead == true) {
                            test->neighbours++;
                        } else {
                            ptr->neighbours++;
                        }
                    } else {
                        matrix_insert(&m, x, b, z, true, 1);
                        ptr = ((struct node*)array->data) + k;
                    }
                }
            }
            i = val_mod(++i, SIZE);
        }

        for (int a1 = 0, i = lower; a1 < matrix_size(id, process, SIZE) + 2; a1++) {
            for (j = 0; j < SIZE; j++) {
                array = get_matrix(&m, i, j);
                if (!array) {
                    continue;
                }
                
                for (int k = (int)array->used - 1; k >= 0; k--) {
                    ptr = ((struct node*)array->data) + k;
                    if (ptr->dead) {
                        if (ptr->neighbours == 2 || ptr->neighbours == 3) {
                            ptr->dead = false;
                        } else {
                            matrix_remove(&m, i, j, ptr->z);
                            ptr = ((struct node*)array->data) + k;
                        }
                    } else {
                        if (ptr->neighbours < 2 || ptr->neighbours > 4) {
                            matrix_remove(&m, i, j, ptr->z);
                            ptr = ((struct node*)array->data) + k;
                        }
                    }
                }
            }
            i = val_mod(++i, SIZE);
        }
    }

    run_time = time_passed + MPI_Wtime() - divide_time;

    if (id == 0) {
        for (int x = 0; x < SIZE; x++) {
            int owner = matrix_owner(x, process, SIZE);
            if (owner != 0) {
                init_recv_row(&m, x, owner);
            }
        }
        matrix_print_live(&m);

        collect_time = time_passed + MPI_Wtime() - divide_time - run_time;

        FILE* out_fp = fopen("time.log", "w");
        char out_str[200];
        sprintf(out_str, "MPI %s: \nscatter_time: %lf \n run_time: %lf\n collect_time: %lf\n   total_time: %lf", input_file, divide_time, run_time, collect_time, (divide_time + run_time + collect_time));
        fwrite(out_str, strlen(out_str), 1, out_fp);
        printf("  execution time: %lf\n", (divide_time + run_time + collect_time));
        } else {
        for (int x = low; x <= high; x++) {
            init_send_row(&m, x, 0);
        }
    }
    MPI_Finalize();
    return 0;
}
