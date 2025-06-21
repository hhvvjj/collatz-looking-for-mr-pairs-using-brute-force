/*
#########################################################################################################
# looking_for_mr_pairs_using_brute_forca.c
#
# Implementation of bruteforce approach to search mr values present on tuple-based transform
#
# The complete article can be found on https://doi.org/10.5281/zenodo.15546925
#
#########################################################################################################
*/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdint.h>
#include <stdbool.h>
#include <omp.h>

/* ============================================================================
 * CONSTANTS AND CONFIGURATION
 * ============================================================================ */

#define MAX_SEQUENCE_LENGTH 50000
#define INITIAL_M_CAPACITY 100
#define HASH_TABLE_SIZE 8192
#define PROGRESS_UPDATE_INTERVAL 3.0
#define PROGRESS_CHECK_FREQUENCY 100000
#define MIN_CHUNK_SIZE 100
#define MAX_CHUNK_SIZE 10000
#define CHUNK_MULTIPLIER 10

/* ============================================================================
 * DATA STRUCTURES
 * ============================================================================ */

// Hash table node for fast m value lookup
typedef struct HashNode {
    uint64_t value;
    struct HashNode* next;
} HashNode;

// Container for m values with hash table optimization
typedef struct {
    HashNode* buckets[HASH_TABLE_SIZE];
    uint64_t* values;
    int count;
    int capacity;
} MValues;

// Thread-safe set for unique mr values
typedef struct {
    uint64_t* values;
    uint64_t* first_n;
    int count;
    int capacity;
    omp_lock_t lock;
} UniqueMrSet;

// Progress tracking structure
typedef struct {
    uint64_t processed;
    uint64_t found_count;
    uint64_t last_n_with_new_unique;
    double last_update_time;
    omp_lock_t lock;
} ProgressTracker;

// Main search context containing all state
typedef struct {
    uint64_t max_n;
    int exponent;
    UniqueMrSet* unique_set;
    ProgressTracker* progress;
    double start_time;
} SearchContext;

// Scheduling configuration
typedef enum {
    SCHEDULE_STATIC = 0,
    SCHEDULE_GUIDED = 1,
    SCHEDULE_DYNAMIC = 2
} SchedulingType;

typedef struct {
    SchedulingType type;
    uint64_t chunk_size;
    int num_threads;
} SchedulingConfig;

/* ============================================================================
 * UTILITY FUNCTIONS
 * ============================================================================ */

// Hash function for fast lookup
static inline uint64_t hash_function(uint64_t value) {
    return value & (HASH_TABLE_SIZE - 1);
}

// Calculate m parameter for a Collatz value
static inline uint64_t calculate_m(uint64_t c) {
    uint64_t p = (c & 1) ? 1 : 2;
    return (c - p) >> 1;
}

// Apply Collatz transformation with overflow check
static inline bool apply_collatz_transform(uint64_t* n) {
    if (*n & 1) {  // Odd number
        // Check for overflow in 3*n + 1
        if (*n > (UINT64_MAX - 1) / 3) {
            return false;
        }
        uint64_t temp = 3 * (*n);
        if (temp > UINT64_MAX - 1) {
            return false;
        }
        *n = temp + 1;
    } else {  // Even number
        *n = *n >> 1;
    }
    return true;
}

// Safe memory allocation with error checking
static void* safe_malloc(size_t size, const char* context) {
    void* ptr = malloc(size);
    if (!ptr) {
        fprintf(stderr, "Error: Memory allocation failed for %s\n", context);
        exit(1);
    }
    return ptr;
}

// Safe memory reallocation with error checking
static void* safe_realloc(void* ptr, size_t size, const char* context) {
    void* new_ptr = realloc(ptr, size);
    if (!new_ptr) {
        fprintf(stderr, "Error: Memory reallocation failed for %s\n", context);
        exit(1);
    }
    return new_ptr;
}

/* ============================================================================
 * PROGRESS TRACKER IMPLEMENTATION
 * ============================================================================ */

static ProgressTracker* create_progress_tracker(void) {
    ProgressTracker* tracker = safe_malloc(sizeof(ProgressTracker), "progress tracker");
    tracker->processed = 0;
    tracker->found_count = 0;
    tracker->last_n_with_new_unique = 0;
    tracker->last_update_time = 0.0;
    omp_init_lock(&tracker->lock);
    return tracker;
}

static void destroy_progress_tracker(ProgressTracker* tracker) {
    if (tracker) {
        omp_destroy_lock(&tracker->lock);
        free(tracker);
    }
}

static void update_last_n_with_new_unique(ProgressTracker* tracker, uint64_t n) {
    omp_set_lock(&tracker->lock);
    if (n > tracker->last_n_with_new_unique) {
        tracker->last_n_with_new_unique = n;
    }
    omp_unset_lock(&tracker->lock);
}

static void update_progress_if_needed(const SearchContext* ctx) {
    ProgressTracker* tracker = ctx->progress;
    omp_set_lock(&tracker->lock);
    
    double current_time = omp_get_wtime();
    
    if (current_time - tracker->last_update_time >= PROGRESS_UPDATE_INTERVAL) {
        tracker->last_update_time = current_time;
        
        double elapsed = current_time - ctx->start_time;
        double rate = tracker->processed / elapsed;
        double eta = (ctx->max_n - tracker->processed) / rate;
        double progress_percent = (double)tracker->processed / ctx->max_n * 100.0;
        
        printf("Progress: (%.1f%%) | Processed: %lu | Unique values of mr found: %d | %.1f nums/sec | ETA: %.1f min\n",
               progress_percent, tracker->processed, ctx->unique_set->count, rate, eta/60.0);
        fflush(stdout);
    }
    
    omp_unset_lock(&tracker->lock);
}

static void increment_progress_counters(ProgressTracker* tracker, bool mr_found) {
    #pragma omp atomic
    tracker->processed++;
    
    if (mr_found) {
        #pragma omp atomic
        tracker->found_count++;
    }
}

/* ============================================================================
 * M VALUES CONTAINER IMPLEMENTATION
 * ============================================================================ */

static void init_m_values(MValues* mv) {
    mv->capacity = INITIAL_M_CAPACITY;
    mv->values = safe_malloc(mv->capacity * sizeof(uint64_t), "m_values array");
    mv->count = 0;
    
    // Initialize hash table
    for (int i = 0; i < HASH_TABLE_SIZE; i++) {
        mv->buckets[i] = NULL;
    }
}

static void destroy_m_values(MValues* mv) {
    if (!mv) return;
    
    // Clean up hash table
    for (int i = 0; i < HASH_TABLE_SIZE; i++) {
        HashNode* current = mv->buckets[i];
        while (current) {
            HashNode* next = current->next;
            free(current);
            current = next;
        }
        mv->buckets[i] = NULL;
    }
    
    // Clean up values array
    free(mv->values);
    mv->values = NULL;
    mv->count = 0;
    mv->capacity = 0;
}

static bool is_m_repeated(const MValues* mv, uint64_t m) {
    uint64_t hash = hash_function(m);
    HashNode* current = mv->buckets[hash];
    
    while (current) {
        if (current->value == m) {
            return true;
        }
        current = current->next;
    }
    return false;
}

static void add_m_value(MValues* mv, uint64_t m) {
    // Expand array if necessary
    if (mv->count >= mv->capacity) {
        int new_capacity = mv->capacity * 2;
        mv->values = safe_realloc(mv->values, new_capacity * sizeof(uint64_t), "m_values expansion");
        mv->capacity = new_capacity;
    }
    
    mv->values[mv->count++] = m;
    
    // Add to hash table
    uint64_t hash = hash_function(m);
    HashNode* new_node = safe_malloc(sizeof(HashNode), "hash node");
    new_node->value = m;
    new_node->next = mv->buckets[hash];
    mv->buckets[hash] = new_node;
}

/* ============================================================================
 * UNIQUE MR SET IMPLEMENTATION
 * ============================================================================ */

static UniqueMrSet* create_unique_mr_set(void) {
    UniqueMrSet* set = safe_malloc(sizeof(UniqueMrSet), "unique mr set");
    set->capacity = 10000;
    set->values = safe_malloc(set->capacity * sizeof(uint64_t), "unique mr values");
    set->first_n = safe_malloc(set->capacity * sizeof(uint64_t), "unique mr first_n");
    set->count = 0;
    omp_init_lock(&set->lock);
    return set;
}

static void destroy_unique_mr_set(UniqueMrSet* set) {
    if (set) {
        free(set->values);
        free(set->first_n);
        omp_destroy_lock(&set->lock);
        free(set);
    }
}

static bool is_mr_already_found(const UniqueMrSet* set, uint64_t mr) {
    omp_set_lock((omp_lock_t*)&set->lock);
    bool found = false;
    for (int i = 0; i < set->count; i++) {
        if (set->values[i] == mr) {
            found = true;
            break;
        }
    }
    omp_unset_lock((omp_lock_t*)&set->lock);
    return found;
}

static bool add_unique_mr(UniqueMrSet* set, uint64_t mr, uint64_t n) {
    omp_set_lock(&set->lock);
    
    // Check if it already exists
    for (int i = 0; i < set->count; i++) {
        if (set->values[i] == mr) {
            omp_unset_lock(&set->lock);
            return false;
        }
    }
    
    // Expand capacity, if necessary
    if (set->count >= set->capacity) {
        int new_capacity = set->capacity * 2;
        set->values = safe_realloc(set->values, new_capacity * sizeof(uint64_t), "unique mr values expansion");
        set->first_n = safe_realloc(set->first_n, new_capacity * sizeof(uint64_t), "unique mr first_n expansion");
        set->capacity = new_capacity;
    }
    
    set->values[set->count] = mr;
    set->first_n[set->count] = n;
    set->count++;
    
    omp_unset_lock(&set->lock);
    return true;
}

static void report_new_unique_mr(uint64_t mr, uint64_t n, const UniqueMrSet* set, ProgressTracker* tracker) {
    update_last_n_with_new_unique(tracker, n);
    
    #pragma omp critical(discovery_report)
    {
        printf(" [*] New unique value of mr = %lu generated by n = %lu (total unique: %d)\n", 
               mr, n, set->count);
        fflush(stdout);
    }
}

/* ============================================================================
 * COLLATZ SEQUENCE ANALYSIS
 * ============================================================================ */

static uint64_t find_first_mr_in_sequence(uint64_t n_start, bool* found) {
    uint64_t n = n_start;
    MValues m_values;
    init_m_values(&m_values);
    
    uint64_t first_mr = 0;
    *found = false;
    
    // Generate Collatz sequence and search for repetitions in m values
    for (int step = 0; step < MAX_SEQUENCE_LENGTH && n != 1; step++) {
        uint64_t m = calculate_m(n);
        
        // Check repetition BEFORE adding
        if (is_m_repeated(&m_values, m)) {
            first_mr = m;
            *found = true;
            break;
        }
        
        add_m_value(&m_values, m);
        
        // Apply Collatz transformation with overflow check
        if (!apply_collatz_transform(&n)) {
            *found = false;
            break;
        }
        
        // Safety check for very long sequences
        if (step > 10000 && n > n_start * 1000) {
            first_mr = 0;
            *found = true;
            break;
        }
    }
    
    // If we naturally reach n=1 (trivial cycle), then mr=0
    if (n == 1 && !*found) {
        first_mr = 0;
        *found = true;
    }
    
    destroy_m_values(&m_values);
    return first_mr;
}

static void process_single_number(uint64_t n, SearchContext* ctx, uint64_t* local_found, uint64_t* local_processed) {
    bool mr_found = false;
    uint64_t mr = find_first_mr_in_sequence(n, &mr_found);
    
    if (mr_found) {
        (*local_found)++;
        increment_progress_counters(ctx->progress, true);
        
        // Add to unique set and check if new
        bool is_new_unique = add_unique_mr(ctx->unique_set, mr, n);
        
        // Report new unique immediately
        if (is_new_unique) {
            report_new_unique_mr(mr, n, ctx->unique_set, ctx->progress);
        }
    } else {
        increment_progress_counters(ctx->progress, false);
    }
    
    (*local_processed)++;
}

/* ============================================================================
 * PARALLEL SCHEDULING CONFIGURATION
 * ============================================================================ */

static SchedulingConfig configure_parallel_scheduling(uint64_t max_n) {
    SchedulingConfig config;
    config.num_threads = omp_get_max_threads();
    
    // Choose scheduling strategy based on problem size
    if (max_n < 1000000) {
        config.type = SCHEDULE_STATIC;
        config.chunk_size = 0; // Let OpenMP decide
    } else if (max_n < 100000000) {
        config.type = SCHEDULE_GUIDED;
        config.chunk_size = 0; // Not applicable for guided
    } else {
        config.type = SCHEDULE_DYNAMIC;
        // Calculate optimal chunk size
        config.chunk_size = max_n / (config.num_threads * CHUNK_MULTIPLIER);
        if (config.chunk_size < MIN_CHUNK_SIZE) config.chunk_size = MIN_CHUNK_SIZE;
        if (config.chunk_size > MAX_CHUNK_SIZE) config.chunk_size = MAX_CHUNK_SIZE;
    }
    
    return config;
}

static void print_scheduling_info(const SchedulingConfig* config) {
    printf("Using scheduling strategy: ");
    switch(config->type) {
        case SCHEDULE_STATIC:
            printf("static (small range)\n");
            break;
        case SCHEDULE_GUIDED:
            printf("guided (medium range)\n");
            break;
        case SCHEDULE_DYNAMIC:
            printf("dynamic with chunk size %lu (large range)\n", config->chunk_size);
            break;
    }
    printf("Number of threads: %d\n\n", config->num_threads);
}

/* ============================================================================
 * MAIN PARALLEL SEARCH EXECUTION
 * ============================================================================ */

static void execute_search_with_static_scheduling(SearchContext* ctx, uint64_t* total_found, uint64_t* total_processed) {
    uint64_t local_found = 0, local_processed = 0;
    
    #pragma omp parallel reduction(+:local_found, local_processed)
    {
        uint64_t thread_found = 0, thread_processed = 0;
        int thread_num = omp_get_thread_num();
        
        #pragma omp for schedule(static)
        for (uint64_t n = 1; n < ctx->max_n; n++) {
            process_single_number(n, ctx, &thread_found, &thread_processed);
            
            if (thread_num == 0 && thread_processed % PROGRESS_CHECK_FREQUENCY == 0) {
                update_progress_if_needed(ctx);
            }
        }
        
        local_found += thread_found;
        local_processed += thread_processed;
    }
    
    *total_found = local_found;
    *total_processed = local_processed;
}

static void execute_search_with_guided_scheduling(SearchContext* ctx, uint64_t* total_found, uint64_t* total_processed) {
    uint64_t local_found = 0, local_processed = 0;
    
    #pragma omp parallel reduction(+:local_found, local_processed)
    {
        uint64_t thread_found = 0, thread_processed = 0;
        int thread_num = omp_get_thread_num();
        
        #pragma omp for schedule(guided)
        for (uint64_t n = 1; n < ctx->max_n; n++) {
            process_single_number(n, ctx, &thread_found, &thread_processed);
            
            if (thread_num == 0 && thread_processed % PROGRESS_CHECK_FREQUENCY == 0) {
                update_progress_if_needed(ctx);
            }
        }
        
        local_found += thread_found;
        local_processed += thread_processed;
    }
    
    *total_found = local_found;
    *total_processed = local_processed;
}

static void execute_search_with_dynamic_scheduling(SearchContext* ctx, uint64_t chunk_size, uint64_t* total_found, uint64_t* total_processed) {
    uint64_t local_found = 0, local_processed = 0;
    
    #pragma omp parallel reduction(+:local_found, local_processed)
    {
        uint64_t thread_found = 0, thread_processed = 0;
        int thread_num = omp_get_thread_num();
        
        #pragma omp for schedule(dynamic, chunk_size)
        for (uint64_t n = 1; n < ctx->max_n; n++) {
            process_single_number(n, ctx, &thread_found, &thread_processed);
            
            if (thread_num == 0 && thread_processed % PROGRESS_CHECK_FREQUENCY == 0) {
                update_progress_if_needed(ctx);
            }
        }
        
        local_found += thread_found;
        local_processed += thread_processed;
    }
    
    *total_found = local_found;
    *total_processed = local_processed;
}

static void execute_parallel_search(SearchContext* ctx, uint64_t* found_count, uint64_t* processed_count) {
    SchedulingConfig config = configure_parallel_scheduling(ctx->max_n);
    print_scheduling_info(&config);
    
    // Execute search with appropriate scheduling
    switch(config.type) {
        case SCHEDULE_STATIC:
            execute_search_with_static_scheduling(ctx, found_count, processed_count);
            break;
        case SCHEDULE_GUIDED:
            execute_search_with_guided_scheduling(ctx, found_count, processed_count);
            break;
        case SCHEDULE_DYNAMIC:
            execute_search_with_dynamic_scheduling(ctx, config.chunk_size, found_count, processed_count);
            break;
    }
}

/* ============================================================================
 * RESULTS OUTPUT AND REPORTING
 * ============================================================================ */

static int* create_sorting_indices(const UniqueMrSet* set) {
    int* indices = safe_malloc(set->count * sizeof(int), "sorting indices");
    
    for (int i = 0; i < set->count; i++) {
        indices[i] = i;
    }
    
    // Simple bubble sort for mr values
    for (int i = 0; i < set->count - 1; i++) {
        for (int j = 0; j < set->count - 1 - i; j++) {
            if (set->values[indices[j]] > set->values[indices[j + 1]]) {
                int temp = indices[j];
                indices[j] = indices[j + 1];
                indices[j + 1] = temp;
            }
        }
    }
    
    return indices;
}

static void write_results_to_file(const SearchContext* ctx, const char* filename) {
    FILE* output = fopen(filename, "w");
    if (!output) {
        printf("Error: Could not create the output file %s\n", filename);
        return;
    }
    
    // Write CSV header and data
    fprintf(output, "n,mr\n");
    for (int i = 0; i < ctx->unique_set->count; i++) {
        fprintf(output, "%lu,%lu\n", ctx->unique_set->first_n[i], ctx->unique_set->values[i]);
    }
    
    fclose(output);
}

static void print_validation_results(const SearchContext* ctx, double total_time, uint64_t found_count, uint64_t processed_count) {
    UniqueMrSet* set = ctx->unique_set;
    int* indices = create_sorting_indices(set);
    
    printf("\nVALIDATION RESULTS:\n");
    printf("\nList of unique mr values found and associated n value:\n\n");
    
    for (int i = 0; i < set->count; i++) {
        int idx = indices[i];
        printf(" [*] mr = %lu (generated using n = %lu)\n",
               set->values[idx], set->first_n[idx]);
    }
    
    printf("\nFinal Summary:\n\n");
    printf(" [*] Exponent: %d\n", ctx->exponent);
    printf(" [*] Total time: %.2f seconds\n", total_time);
    printf(" [*] Speed: %.1f numbers/second\n", (double)processed_count / total_time);
    printf(" [*] Numbers processed: 2^%d = %lu\n", ctx->exponent, processed_count);
    printf(" [*] Numbers with mr found: %lu\n", found_count);
    printf(" [*] Percentage with mr: %.2f%%\n", (double)found_count / processed_count * 100.0);
    printf(" [*] Unique mr values found: %lu\n", (uint64_t)set->count);
    
    // Show sorted mr values for easy copying
    printf("\nFinal list of mr values:\n\n");
    for (int i = 0; i < set->count; i++) {
        int idx = indices[i];
        if (i > 0) printf(", ");
        printf("%lu", set->values[idx]);
    }
    printf("\n");
    
    free(indices);
}

/* ============================================================================
 * COMMAND LINE ARGUMENT PROCESSING
 * ============================================================================ */

static bool validate_and_parse_arguments(int argc, char* argv[], int* exponent, uint64_t* max_n) {
    if (argc != 2) {
        printf("Usage: %s <exponent>\n", argv[0]);
        printf("Example: %s 25  (to search n < 2^25)\n", argv[0]);
        printf("Recommended exponents:\n");
        printf("  20 -> 2^20 = 1,048,576 (quick test)\n");
        printf("  25 -> 2^25 = 33,554,432 (default)\n");
        printf("  30 -> 2^30 = 1,073,741,824 (intensive use)\n");
        return false;
    }
    
    *exponent = atoi(argv[1]);
    
    if (*exponent < 1 || *exponent > 64) {
        printf("Error: Exponent must be between 1 and 64\n");
        printf("Exponent %d is out of valid range\n", *exponent);
        return false;
    }
    
    *max_n = 1UL << *exponent;
    return true;
}

static void print_program_header(int exponent, uint64_t max_n) {
    printf("\n-----------------------------------------------------------------------------------------\n");
    printf("- Bruteforce Collatz sequences using Tuple-based transform to compute first repetitions -\n");
    printf("-----------------------------------------------------------------------------------------\n");
    printf("\nUsing %d threads for parallelization\n", omp_get_max_threads());
    printf("Searching for the first pair mr between sequences of n < 2^%d...\n", exponent);
    printf("Range: from 1 to %lu \n\n", max_n - 1);
}

/* ============================================================================
 * MAIN FUNCTION
 * ============================================================================ */

int main(int argc, char* argv[]) {
    // Parse and validate command line arguments
    int exponent;
    uint64_t max_n;
    if (!validate_and_parse_arguments(argc, argv, &exponent, &max_n)) {
        return 1;
    }
    
    // Print program header
    print_program_header(exponent, max_n);
    
    // Initialize search context
    SearchContext ctx = {
        .max_n = max_n,
        .exponent = exponent,
        .unique_set = create_unique_mr_set(),
        .progress = create_progress_tracker(),
        .start_time = omp_get_wtime()
    };
    
    // Execute the parallel search
    uint64_t found_count = 0;
    uint64_t processed_count = 0;
    
    execute_parallel_search(&ctx, &found_count, &processed_count);
    
    double end_time = omp_get_wtime();
    double total_time = end_time - ctx.start_time;
    
    // Generate output filename and write results
    char filename[256];
    snprintf(filename, sizeof(filename), "mr_results_2_%d.txt", exponent);
    write_results_to_file(&ctx, filename);
    
    // Print validation results and summary
    print_validation_results(&ctx, total_time, found_count, processed_count);
    
    printf("\n\nResults saved to '%s'\n", filename);
    
    // Cleanup
    destroy_unique_mr_set(ctx.unique_set);
    destroy_progress_tracker(ctx.progress);
    
    return 0;
}