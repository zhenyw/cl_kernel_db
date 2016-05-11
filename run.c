
#define _GNU_SOURCE

#include <fcntl.h>
#include <unistd.h>
#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <stdbool.h>
#include <string.h>
#include <math.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <assert.h>
#include <CL/opencl.h>

#include <sys/time.h>

char *KernelSource = NULL;
char *KernelOption = NULL;

bool get_kernel_src(char *src)
{
    FILE *f;
    long s;
    struct stat st;

    if (stat(src, &st) < 0)
	return false;

    if (st.st_size == 0)
	return false;

    f = fopen(src, "r");
    if (!f)
	return false;
    fseek(f, 0, SEEK_END);
    s = ftell(f);
    fseek(f, 0, SEEK_SET);
    KernelSource = (char *)malloc(sizeof(char)*(s + 1));
    if (!KernelSource) {
	fclose(f);
	return false;
    }
    memset(KernelSource, 0, sizeof(char)*(s + 1));
    fread(KernelSource, sizeof(char), s, f);
    fclose(f);
    return true;
}

bool get_kernel_option(char *opt)
{
    FILE *f;
    long s;
    struct stat st;

    if (stat(opt, &st) < 0)
	return false;

    if (st.st_size == 0)
	return false;

    f = fopen(opt, "r");
    if (!f)
	return false;
    fseek(f, 0, SEEK_END);
    s = ftell(f);
    fseek(f, 0, SEEK_SET);
    KernelOption = (char *)malloc(sizeof(char)*(s+1));
    if (!KernelOption) {
	fclose(f);
	return false;
    }
    memset(KernelOption, 0, sizeof(char)*(s+1));
    fread(KernelOption, sizeof(char), s, f);
    fclose(f);
    return true;
}

/* return millisecond */
double get_time_diff(struct timeval *s, struct timeval *e)
{
    double s_ms , e_ms , diff;

    s_ms = (double)s->tv_sec*1000 + (double)s->tv_usec / 1000.0;
    e_ms = (double)e->tv_sec*1000 + (double)e->tv_usec / 1000.0;

    diff = (double)e_ms - (double)s_ms;

    return diff;
}

static char *
beignet_find_magic(char *base, unsigned int b_length, unsigned int *left,
		   const char *magic, int m_length)
{
    int i, j;

    for (i = 0; i < b_length; i = i + j + 1) {
	for (j = 0; j < m_length; j++)
	    if (base[i+j] != magic[j])
		break;
	if (j == m_length) {
	    *left = b_length - i;
	    return &base[i];
	}
    }
    return NULL;
}

#define BEIGNET_BIN_HEADER_LEN 8

#define KERNEL_MAGIC_SIZE 4
#define KERNEL_NAME_SIZE sizeof(size_t)
#define KERNEL_ARG_NUM_SIZE 4
#define KERNEL_ARG_TYPE_SIZE 4
#define KERNEL_ARG_SIZE_SIZE 4
#define KERNEL_ARG_ALIGN_SIZE 4
#define KERNEL_ARG_BTI_SIZE 1
#define KERNEL_ARG_ADDR_SPACE_SIZE 4
#define KERNEL_ARG_TYPE_NAME_SIZE sizeof(size_t)
#define KERNEL_ARG_ACCESS_QUAL_SIZE sizeof(size_t)
#define KERNEL_ARG_TYPE_QUAL_SIZE sizeof(size_t)
#define KERNEL_ARG_NAME_SIZE sizeof(size_t)
#define KERNEL_PATCHES_SIZE sizeof(size_t)
#define KERNEL_PATCHES_TYPE_SIZE sizeof(unsigned int)
#define KERNEL_PATCHES_SUBTYPE_SIZE sizeof(unsigned int)
#define KERNEL_PATCHES_OFFSET_SIZE sizeof(unsigned int)
#define KERNEL_CURBE_SIZE 4
#define KERNEL_SIMD_SIZE 4
#define KERNEL_STACK_SIZE 4
#define KERNEL_STRATCH_SIZE 4
#define KERNEL_USE_SLM_SIZE 1
#define KERNEL_SLM_SIZE 4
#define KERNEL_COMPILE_WG_SIZE (sizeof(size_t)*3)
#define KERNEL_SAMPLE_SET_SIZE sizeof(int)
#define KERNEL_SAMPLE_MAGIC_SIZE 4
#define KERNEL_SAMPLE_MAP_SIZE sizeof(size_t)
#define KERNEL_SAMPLE_MAP_ITEM_SIZE (4*2)
#define KERNEL_SAMPLE_MAGIC_END_SIZE 4
#define KERNEL_SAMPLE_RET_SIZE sizeof(size_t)
#define KERNEL_IMAGE_SET_SIZE sizeof(int)
#define KERNEL_IMAGE_MAGIC_SIZE 4
#define KERNEL_IMAGE_REGMAP_SIZE sizeof(size_t)
#define KERNEL_IMAGE_MAP_FIRST_SIZE 4 /* Register */
#define KERNEL_IMAGE_MAP_SECOND_SIZE (4*8)
#define KERNEL_IMAGE_IDXMAP_SIZE sizeof(size_t)
#define KERNEL_IMAGE_IDXMAP_FIRST_SIZE 4
#define KERNEL_IMAGE_IDXMAP_SECOND_SIZE (4*8)	    
#define KERNEL_IMAGE_MAGIC_END_SIZE 4
#define KERNEL_IMAGE_RET_SIZE sizeof(size_t)
#define KERNEL_CODE_SIZE sizeof(size_t)
	    
static size_t
beignet_kernel_code_size(char *kern)
{
    uint32_t f_32, i;
    size_t f_size_t;
    int f_int;

    if (kern[0] != 'N' || kern[1] != 'R' || kern[2] != 'E' || kern[3] != 'K') {
	fprintf(stderr, "Wrong beignet kernel section\n");
	return 0;
    }
    
    kern += KERNEL_MAGIC_SIZE;
    f_size_t = *(size_t *)kern;
    kern += KERNEL_NAME_SIZE;
    kern += f_size_t;
    
    f_32 = *(uint32_t *)kern;
    kern += KERNEL_ARG_NUM_SIZE;
    for (i = 0; i < f_32; i++) {
	uint32_t tmp;
	kern += KERNEL_ARG_TYPE_SIZE;
	kern += KERNEL_ARG_SIZE_SIZE;
	kern += KERNEL_ARG_ALIGN_SIZE;
	kern += KERNEL_ARG_BTI_SIZE;
	kern += KERNEL_ARG_ADDR_SPACE_SIZE;

	tmp = *(uint32_t *)kern;
	kern += KERNEL_ARG_TYPE_NAME_SIZE;
	kern += tmp;

	tmp = *(uint32_t *)kern;
	kern += KERNEL_ARG_ACCESS_QUAL_SIZE;
	kern += tmp;

	tmp = *(uint32_t *)kern;
	kern += KERNEL_ARG_TYPE_QUAL_SIZE;
	kern += tmp;

	tmp = *(uint32_t *)kern;
	kern += KERNEL_ARG_NAME_SIZE;
	kern += tmp;
    }

    f_size_t = *(size_t *)kern;
    kern += KERNEL_PATCHES_SIZE;
    for (i = 0; i < f_size_t; i++) {
	kern += (KERNEL_PATCHES_TYPE_SIZE +
		 KERNEL_PATCHES_SUBTYPE_SIZE +
		 KERNEL_PATCHES_OFFSET_SIZE);
    }

    kern += KERNEL_CURBE_SIZE;
    kern += KERNEL_SIMD_SIZE;
    kern += KERNEL_STACK_SIZE;
    kern += KERNEL_STRATCH_SIZE;
    kern += KERNEL_USE_SLM_SIZE;
    kern += KERNEL_SLM_SIZE;
    kern += KERNEL_COMPILE_WG_SIZE;

    f_int = *(int *)kern;
    kern += KERNEL_SAMPLE_SET_SIZE;
    if (f_int != 0) {
	kern += KERNEL_SAMPLE_MAGIC_SIZE;
	f_size_t = *(size_t *)kern;
	kern += KERNEL_SAMPLE_MAP_SIZE;
	for (i = 0; i < f_size_t; i++) {
	    kern += KERNEL_SAMPLE_MAP_ITEM_SIZE;
	}
	kern += KERNEL_SAMPLE_MAGIC_END_SIZE;
	kern += KERNEL_SAMPLE_RET_SIZE;
    }

    f_int = *(int *)kern;
    kern += KERNEL_IMAGE_SET_SIZE;
    if (f_int != 0) {
	kern += KERNEL_IMAGE_MAGIC_SIZE;
	f_size_t = *(size_t *)kern;
	kern += KERNEL_IMAGE_REGMAP_SIZE;
	for (i = 0; i < f_size_t; i++) {
	    kern += KERNEL_IMAGE_MAP_FIRST_SIZE;
	    kern += KERNEL_IMAGE_MAP_SECOND_SIZE;
	}
	f_size_t = *(size_t *)kern;
	kern += KERNEL_IMAGE_IDXMAP_SIZE;
	for (i = 0; i < f_size_t; i++) {
	    kern += KERNEL_IMAGE_IDXMAP_FIRST_SIZE;
	    kern += KERNEL_IMAGE_IDXMAP_SECOND_SIZE;
	}
	kern += KERNEL_IMAGE_MAGIC_END_SIZE;
	kern += KERNEL_IMAGE_RET_SIZE;
    }

    return *(size_t *)kern;
}

void usage()
{
    printf("./run [-p 0/1] <cl kernel file> [<cl kernel option file>]\n");
    printf("\n");
    printf("-p : 0 uses beignet platform, 1 uses VPG driver platform.\n");
    exit(0);
}

int main(int argc, char *argv[])
{
    cl_int err;
    cl_platform_id platform_ids[2];
    cl_uint num_platform;
    bool found_platform = false;
    cl_platform_id target_platform_id;
    cl_device_id device_id;
    cl_context context;
    cl_command_queue commands;
    cl_program program;
    char platname[64];
    char drvname[64];
    size_t ret;
    struct timeval start, end;
    char *cl_src, *cl_opt = NULL;
    int platform_select = 0;            /* 0: beignet, 1: VPG */
    int opt_idx, i;

    if (argc < 2)
	usage();

    opt_idx = 1;
    if (strcmp(argv[opt_idx], "-p") == 0) {
	if (strcmp(argv[opt_idx+1], "0") == 0)
	    platform_select = 0;
	else if (strcmp(argv[opt_idx+1], "1") == 0)
	    platform_select = 1;
	else
	    usage();
	opt_idx += 2;
    }

    cl_src = argv[opt_idx];
    if (argv[opt_idx+1])
	cl_opt = argv[opt_idx+1];
    
    err = clGetPlatformIDs(0, NULL, &num_platform);
    assert(err == CL_SUCCESS);

    if (num_platform > 2) {
	fprintf(stderr, "Invalid platform\n");
	exit(1);
    }

    err = clGetPlatformIDs(num_platform, platform_ids, NULL);
    assert(err == CL_SUCCESS);

    for (i = 0; i < num_platform; i++) {
	memset(platname, 0, sizeof(platname));
	err = clGetPlatformInfo(platform_ids[i], CL_PLATFORM_NAME, 64, platname, &ret);
	assert(err == CL_SUCCESS);
	platname[ret] = '\0';

	if ((platform_select == 0 && strcmp(platname, "Intel Gen OCL Driver") == 0) ||
	    (platform_select == 1 && strcmp(platname, "Intel(R) OpenCL") == 0)) {
	    target_platform_id = platform_ids[i];
	    found_platform = true;
	    break;
	}
    }

    if (!found_platform) {
	fprintf(stderr, "No valid platform found\n");
	exit(1);
    }

    err = clGetDeviceIDs(target_platform_id, CL_DEVICE_TYPE_GPU, 1, &device_id, NULL);
    assert(err == CL_SUCCESS);

    context = clCreateContext(NULL, 1, &device_id, NULL, NULL, &err);
    assert(context && err == CL_SUCCESS);

#if 0    
    err = clGetDeviceInfo(device_id, CL_DRIVER_VERSION, 128, drvname, &ret);
    assert(err == CL_SUCCESS);
    drvname[ret] = '\0';
    printf("Driver: %s\n", drvname);
#endif
    
    commands = clCreateCommandQueue(context, device_id, 0, &err);
    assert(commands && err == CL_SUCCESS);

    if (get_kernel_src(cl_src) == false) {
	fprintf(stderr, "Failed to open kernel file\n");
	exit(1);
    }

    if (cl_opt && get_kernel_option(cl_opt) == false) {
	fprintf(stderr, "Invalid kernel option file\n");
	exit(1);
    }

    program = clCreateProgramWithSource(context, 1, (const char **) & KernelSource, NULL, &err);
    assert(program && err == CL_SUCCESS);

    gettimeofday(&start, NULL);
	
    err = clBuildProgram(program, 1, &device_id, KernelOption ? KernelOption : NULL, NULL, NULL);
    if (err != CL_SUCCESS)
    {
        size_t len;
        char buffer[2048];

        fprintf(stderr, "Error: Failed to build program executable!\n");
        clGetProgramBuildInfo(program, device_id, CL_PROGRAM_BUILD_LOG, sizeof(buffer), buffer, &len);
        fprintf(stderr, "%s\n", buffer);
        exit(1);
    }
    gettimeofday(&end, NULL);

    size_t bin_size;
    err = clGetProgramInfo(program, CL_PROGRAM_BINARY_SIZES, sizeof(bin_size), &bin_size, NULL);
    assert(err == CL_SUCCESS);

    char *bin = (char *)malloc(bin_size);
    if (!bin) {
	fprintf(stderr, "Error: binary buffer alloc fail\n");
	exit(1);
    }
    err = clGetProgramInfo(program, CL_PROGRAM_BINARIES, bin_size, &bin, NULL);
    assert(err == CL_SUCCESS);

#if 0
    FILE *f_bin = fopen("kernel.bin", "w");
    fwrite(bin, bin_size, 1, f_bin);
    fflush(f_bin);
    fclose(f_bin);
#endif
    
    unsigned int code_size = 0;
    
    if (platform_select == 0) {
	char *new_kernel = bin;
	unsigned int len = bin_size;
	unsigned int left;

	if (bin[0] != 0 && bin[1] != 'G' && bin[2] != 'E' &&
	    bin[3] != 'N' && bin[4] != 'C') {
	    fprintf(stderr, "Error: invalid beignet binary format\n");
	    exit(1);
	}

	new_kernel += BEIGNET_BIN_HEADER_LEN;
	
	while ((new_kernel = beignet_find_magic(new_kernel, len, &left, "NREK", 4)) != NULL) {
	    len = left;
	    code_size += beignet_kernel_code_size(new_kernel);
	    new_kernel += 4;
	}
    }

    printf("%s,%s,%u,%u,%.2f\n", cl_src, platname, bin_size, code_size, get_time_diff(&start, &end));
    
    clReleaseProgram(program);
    clReleaseCommandQueue(commands);
    clReleaseContext(context);

    return 0;
}
