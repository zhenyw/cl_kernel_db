__kernel void foo(__global unsigned char* inBuf, __global unsigned char* outBuf)
{
    size_t idx = get_global_id(0);
    outBuf[idx] = inBuf[idx] + 1;
}
