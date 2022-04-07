#include <iostream> 
#include <cassert>
#include <queue>
#include "../include/intersection.h"

using namespace std;

extern bool calculate_domains_context;

#ifdef PRINT_INTERSECTION
size_t CD_totalEdgeCase1 = 0;
size_t CD_totalEdgeCase2 = 0;
size_t CD_totalLenRare = 0;
size_t CD_totalLenFreq = 0;
size_t CD_totalGalloping = 0;
size_t CD_totalA = 0;
size_t CD_totalB = 0;
size_t CD_totalShuffling = 0;
size_t CD_totalCountGalloping = 0;
size_t CD_totalCountShuffling = 0;

size_t totalEdgeCase1 = 0;
size_t totalEdgeCase2 = 0;
size_t totalLenRare = 0;
size_t totalLenFreq = 0;
size_t totalGalloping = 0;
size_t totalA = 0;
size_t totalB = 0;
size_t totalShuffling = 0;
size_t totalCountGalloping = 0;
size_t totalCountShuffling = 0;
#endif

#ifdef COMPRESS

void quit()
{
    system("pause");
    exit(0);
}

void align_malloc(void **memptr, size_t alignment, size_t size)
{
    int malloc_flag = posix_memalign(memptr, alignment, size);
    if (malloc_flag) {
        std::cerr << "posix_memalign: " << strerror(malloc_flag) << std::endl;
        quit();
    }
}

uint32_t* prepare_shuffling_dict_avx()
{
    uint32_t* arr = new uint32_t[2048];
    for(int i=0; i<256; ++i){
        int count=0, rest=7;
        for(int b=0; b<8; ++b){
            if(i & (1 << b)){
                // n index at pos p - move nth element to pos p
                arr[i*8 + count] = b; // move all set bits to beginning
                ++count;
            }else{
                arr[i*8 + rest] = b; // move rest at the end
                --rest;
            }
        }
    }
    return arr;
}

int* prepare_byte_check_mask_dict2()
{
    int * mask = new int[65536];

    auto trans_c_s = [](const int c) -> int {
        switch (c) {
            case 0: return -1; // no match
            case 1: return 0;
            case 2: return 1;
            case 4: return 2;
            case 8: return 3;
            default: return 4; // multiple matches.
        }
    };

    for (int x = 0; x < 65536; ++x) {        
        int c0 = (x & 0xf), c1 = ((x >> 4) & 0xf);
        int c2 = ((x >> 8) & 0xf), c3 = ((x >> 12) & 0xf);
        int s0 = trans_c_s(c0), s1= trans_c_s(c1);
        int s2 = trans_c_s(c2), s3 = trans_c_s(c3);
        
        bool is_multiple_match = (s0 == 4) || (s1 == 4) ||
                (s2 == 4) || (s3 == 4);
        if (is_multiple_match) {
            mask[x] = -1;
            continue;
        }
        bool is_no_match = (s0 == -1) && (s1 == -1) &&
                (s2 == -1) && (s3 == -1);
        if (is_no_match) {
            mask[x] = -2;
            continue;
        }
        if (s0 == -1) s0 = 0; if (s1 == -1) s1 = 1;
        if (s2 == -1) s2 = 2; if (s3 == -1) s3 = 3;
        mask[x] = (s0) | (s1 << 2) | (s2 << 4) | (s3 << 6);        
    }

    return mask;
}

uint8_t * prepare_match_shuffle_dict2()
{
    uint8_t * dict = new uint8_t[4096];

    for (int x = 0; x < 256; ++x) {
        for (int i = 0; i < 4; ++i) {
            uint8_t c = (x >> (i << 1)) & 3; // c = 0, 1, 2, 3
            int pos = x * 16 + i * 4;
            for (uint8_t j = 0; j < 4; ++j)
                dict[pos + j] = c * 4 + j;
        }
    }

    return dict;
}

uint8_t* prepare_shuffling_dict_u16()
{
    uint8_t* mask = new uint8_t[4096];
    memset(mask, 0xff, sizeof(uint8_t) * 4096);
    int size = 0;
    for (int i = 0; i < 256; ++i) {
        int counter = 0;
        for (int j = 0; j < 8; ++j) {
            if (i & (1 << j)) {
                mask[size + counter * 2    ] = 2 * j;
                mask[size + counter * 2 + 1] = 2 * j + 1;
                ++counter;             
             }              
        }
        size += 16;
    }
    return mask;
}

int intersect_qfilter_bsr_b4_v2(const int* bases_a, const int* states_a, int size_a,
            const int* bases_b, const int* states_b, int size_b,
            int* bases_c, int* states_c) //, int& num_int_c)
{
    int i = 0, j = 0, size_c = 0; //, int_c = 0;
    int qs_a = size_a - (size_a & 3);
    int qs_b = size_b - (size_b & 3);

    __m128i base_a = _mm_lddqu_si128((__m128i*)bases_a + i);
    __m128i base_b = _mm_lddqu_si128((__m128i*)bases_b + j);
    __m128i byte_group_a = _mm_shuffle_epi8(base_a, byte_check_group_a_order[0]);
    __m128i byte_group_b = _mm_shuffle_epi8(base_b, byte_check_group_b_order[0]);
    __m128i cmp_mask, and_state;

    //std::cout << "load bases done" << std::endl;

    while (i < qs_a && j < qs_b) {
        __m128i byte_check_mask = _mm_cmpeq_epi8(byte_group_a, byte_group_b);
        int bc_mask = _mm_movemask_epi8(byte_check_mask);
        int ms_order = byte_check_mask_dict[bc_mask];
        //std::cout << "order = " << ms_order << std::endl;

        if (__builtin_expect(ms_order != -2, 0)) {
            //std::cout << "state a = " << states_a[i] << std::endl;
            //std::cout << "state b = " << states_b[j] << std::endl;
            __m128i state_a = _mm_lddqu_si128((__m128i*)(states_a + i));
            __m128i state_b = _mm_lddqu_si128((__m128i*)(states_b + j));    
            //std::cout << "load states done" << std::endl;
            if (ms_order > 0) {                
                __m128i sf_base_b = _mm_shuffle_epi8(base_b, match_shuffle_dict[ms_order]);
                __m128i sf_state_b = _mm_shuffle_epi8(state_b, match_shuffle_dict[ms_order]);
                cmp_mask = _mm_cmpeq_epi32(base_a, sf_base_b);                
                and_state = _mm_and_si128(state_a, sf_state_b);
                __m128i state_mask = _mm_cmpeq_epi32(and_state, all_zero_si128);
                cmp_mask = _mm_andnot_si128(state_mask, cmp_mask);
            } else {
                __m128i cmp_mask0 = _mm_cmpeq_epi32(base_a, base_b);
                __m128i state_c0 = _mm_and_si128(
                        _mm_and_si128(state_a, state_b), cmp_mask0);
                __m128i base_sf1 = _mm_shuffle_epi32(base_b, cyclic_shift1);
                __m128i state_sf1 = _mm_shuffle_epi32(state_b, cyclic_shift1);
                __m128i cmp_mask1 = _mm_cmpeq_epi32(base_a, base_sf1);
                __m128i state_c1 = _mm_and_si128(
                        _mm_and_si128(state_a, state_sf1), cmp_mask1);
                __m128i base_sf2 = _mm_shuffle_epi32(base_b, cyclic_shift2);
                __m128i state_sf2 = _mm_shuffle_epi32(state_b, cyclic_shift2);
                __m128i cmp_mask2 = _mm_cmpeq_epi32(base_a, base_sf2);
                __m128i state_c2 = _mm_and_si128(
                        _mm_and_si128(state_a, state_sf2), cmp_mask2);
                __m128i base_sf3 = _mm_shuffle_epi32(base_b, cyclic_shift3);
                __m128i state_sf3 = _mm_shuffle_epi32(state_b, cyclic_shift3);
                __m128i cmp_mask3 = _mm_cmpeq_epi32(base_a, base_sf3);
                __m128i state_c3 = _mm_and_si128(
                        _mm_and_si128(state_a, state_sf3), cmp_mask3);
                and_state = _mm_or_si128(
                        _mm_or_si128(state_c0, state_c1),
                        _mm_or_si128(state_c2, state_c3)
                        );                
                __m128i state_mask = _mm_cmpeq_epi32(and_state, all_zero_si128);
                cmp_mask = _mm_andnot_si128(state_mask, all_one_si128);
            } 

            int mask = _mm_movemask_ps((__m128)cmp_mask);
            __m128i res_b = _mm_shuffle_epi8(base_a, shuffle_mask[mask]);
            __m128i res_s = _mm_shuffle_epi8(and_state, shuffle_mask[mask]);
            _mm_storeu_si128((__m128i*)(bases_c + size_c), res_b);
            _mm_storeu_si128((__m128i*)(states_c + size_c), res_s);

            size_c += _mm_popcnt_u32(mask);        
        }
        int a_max = bases_a[i + 3];
        int b_max = bases_b[j + 3];
        if (a_max <= b_max) {
            i += 4;
            //base_a = _mm_load_si128((__m128i*)(bases_a + i));
            base_a = _mm_loadu_si128((__m128i*)(bases_a + i));
            byte_group_a = _mm_shuffle_epi8(base_a, byte_check_group_a_order[0]);
            _mm_prefetch((char*) (bases_a + i + 16), _MM_HINT_T0);
            _mm_prefetch((char*) (states_a + i + 16), _MM_HINT_T0);
        }
        if (a_max >= b_max) {
            j += 4;
            //base_b = _mm_load_si128((__m128i*)(bases_b + j));
            base_b = _mm_loadu_si128((__m128i*)(bases_b + j));
            byte_group_b = _mm_shuffle_epi8(base_b, byte_check_group_b_order[0]);
            _mm_prefetch((char*) (bases_b + j + 16), _MM_HINT_T0);
            _mm_prefetch((char*) (states_b + j + 16), _MM_HINT_T0);
        }
    }

    while (i < size_a && j < size_b) {
        if (bases_a[i] == bases_b[j]) {
            states_c[size_c] = states_a[i] & states_b[j];
            if (states_c[size_c] != 0) bases_c[size_c++] = bases_a[i];
            i++; j++;
        } else if (bases_a[i] < bases_b[j]){
            i++;
        } else {
            j++;
        }
    }
#ifdef PRINT_INTERSECTION
    if (calculate_domains_context) {
        CD_totalA += size_a;
        CD_totalB += size_b;
        CD_totalCountShuffling += size_c;
        CD_totalShuffling++;
    }
    else {
        totalA += size_a;
        totalB += size_b;
        totalCountShuffling += size_c;
        totalShuffling++;
    }
#endif

    return size_c;
}

int intersect_simdgalloping_bsr(const int* bases_a, const int* states_a, int size_a,
        const int* bases_b, const int* states_b, int size_b,
        int* bases_c, int* states_c)
{
    int i = 0, j = 0, size_c = 0;
    int qs_b = size_b - (size_b & 3);
    for (i = 0; i < size_a; ++i) {
        // double-jump:
        int r = 1;
        while (j + (r << 2) < qs_b && bases_a[i] > bases_b[j + (r << 2) + 3]) r <<= 1;
        // binary search:
        int upper = (j + (r << 2) < qs_b) ? (r) : ((qs_b - j - 4) >> 2);
        if (bases_b[j + (upper << 2) + 3] < bases_a[i]) break;        
        int lower = (r >> 1);
        while (lower < upper) {
            int mid = (lower + upper) >> 1; //10 ms
            if (bases_b[j + (mid << 2) + 3] >= bases_a[i]) upper = mid; //30 ms
            else lower = mid + 1;
        }
        j += (lower << 2);

        __m128i bv_a = _mm_set_epi32(bases_a[i], bases_a[i], bases_a[i], bases_a[i]);
        __m128i bv_b = _mm_lddqu_si128((__m128i*)(bases_b + j));
        __m128i cmp_mask = _mm_cmpeq_epi32(bv_a, bv_b);
        int mask = _mm_movemask_ps((__m128)cmp_mask);
        if (mask != 0) {
            int p = __builtin_ctz(mask);
            states_c[size_c] = states_a[i] & states_b[j + p]; //20 ms
            if (states_c[size_c] != 0) bases_c[size_c++] = bases_a[i];
        }
    }

    while (i < size_a && j < size_b) {
        if (bases_a[i] == bases_b[j]) {
            states_c[size_c] = states_a[i] & states_b[j];
            if (states_c[size_c] != 0) bases_c[size_c++] = bases_a[i];
            i++; j++;
        } else if (bases_a[i] < bases_b[j]){
            i++;
        } else {
            j++;
        }
    }

#ifdef PRINT_INTERSECTION
    if (calculate_domains_context) {
        CD_totalLenRare += size_a;
        CD_totalLenFreq += size_b;
        CD_totalCountGalloping += size_c;
        CD_totalGalloping++;
    }
    else {
        totalLenRare += size_a;
        totalLenFreq += size_b;
        totalCountGalloping += size_c;
        totalGalloping++;
    }
#endif

#ifdef CALCULATE_COMPRESSION
#endif

    return size_c;
}

int offline_uint_trans_bsr(int *set_a, int size_a, int *bases_a, int *states_a)
{
    int cnt = -1;
    for (int i = 0; i < size_a; ++i) {
        int u = set_a[i];
        int u_base = (u >> BSR_SHIFT);
        int u_bit = (1 << (u & BSR_MASK));
        if (cnt == -1 || bases_a[cnt] != u_base) {
            bases_a[++cnt] = u_base;
            states_a[cnt] = u_bit;
        } else {
            states_a[cnt] |= u_bit;
        }
    }
    return ++cnt;    
}

int offline_bsr_trans_uint(int *bases_a, int *states_a, int size_a, int *set_a)
{
    int cnt = 0;
    for (int i = 0; i < size_a; ++i) {
        int u_high = (bases_a[i] << BSR_SHIFT);
        int state = states_a[i];
        while (state) {
            int u = (u_high | __builtin_ctz(state));
            set_a[cnt++] = u;
            state &= (state - 1);
        }
    }
    return cnt;
}



///////////////////////////////////////// APIs ////////////////////////////////////////


//assume size_a <= size_b 
int intersect(const int* bases_a, const int* states_a, int size_a,
              const int* bases_b, const int* states_b, int size_b,
              int* bases_c, int* states_c)
{
	if(size_a == 0 || size_b == 0)
    {
        return 0;
    }
	else if(bases_a[size_a-1] < bases_b[0] || bases_b[size_b-1] < bases_a[0])
    {
        return 0;
    }
    else 
    {
        if(32*size_a < size_b)
        {
            return intersect_simdgalloping_bsr(bases_a, states_a, size_a,
                                               bases_b, states_b, size_b, 
                                               bases_c, states_c);
        }
        else
        {
            return intersect_qfilter_bsr_b4_v2(bases_a, states_a, size_a,
                                               bases_b, states_b, size_b, 
                                               bases_c, states_c);
        }
    }
}

#else

int getBit(int value, int position) {
	return ( ( value & (1 << position) ) >> position);
}

void prepare_shuffling_dictionary() {
	for(int i = 0; i < 16; i++) {
		int counter = 0;
		char permutation[16];
		memset(permutation, 0xFF, sizeof(permutation));
		for(char b = 0; b < 4; b++) {
			if(getBit(i, b)) {
				permutation[counter++] = 4*b;
				permutation[counter++] = 4*b + 1;
				permutation[counter++] = 4*b + 2;
				permutation[counter++] = 4*b + 3;
			}
		}
		__m128i mask = _mm_loadu_si128((const __m128i*)permutation);
		shuffle_mask[i] = mask;
	}
}

/**
 * fast scalar scheme designed by n. kurz.
 */
//s_a > 0 and s_b > 0
inline size_t scalar(const int32_t *a, const size_t lena,
		const int32_t *b, const size_t lenb, int32_t *out)
{	
	const int32_t *starta = a;
	const int32_t *startb = b;
	const int32_t *enda = a + lena;
	const int32_t *endb = b + lenb;
	size_t count = 0;

	if(lena == 0 || lenb == 0)
		return 0;

	while (1) {
		while (*a < *b) {
skip_first_compare:
			if (++a == enda)
				return count;
		}
		while (*a > *b) {
			if (++b == endb)
				return count;
		}
		if (*a == *b) {
			out[count++] = *a;
			if (++a == enda || ++b == endb)
				return count;
		} else {
			goto skip_first_compare;
		}
	}

	return count;
}


int intersect_galloping(const int32_t* rare, const int32_t* freq, size_t lenRare, size_t lenFreq, int* C)
{
#ifdef PRINT_INTERSECTION
    if (calculate_domains_context) {
        CD_totalLenRare += lenRare;
        CD_totalLenFreq += lenFreq;
        CD_totalGalloping++;
    }
    else {
        totalLenRare += lenRare;
        totalLenFreq += lenFreq;
        totalGalloping++;
    }
#endif
	
    size_t count = 0;

#ifndef SIMD_GALLOPING
    //assert(false);
    int first = 0;
    for(size_t i = 0; i < lenRare; ++i)
    {
        int last = lenFreq - 1;
        while(first <= last)
        {
            int mid = (first + last) / 2;
            if(freq[mid] > rare[i])
                last = mid - 1;
            else
                first = mid + 1;
        }
        if(last >= 0 && freq[last] == rare[i]) {
            C[count++] = rare[i];
        }
    }
#ifdef PRINT_INTERSECTION
    if (calculate_domains_context)
        CD_totalCountGalloping += count;
    else
        totalCountGalloping += count;
#endif
    return count;
#else
	int32_t *out = C;

	const int32_t * const startRare = rare;
	const int32_t * const startFreq = freq;

	typedef __m128i vec;
	const int32_t veclen = sizeof(vec) / sizeof(int32_t);
	const size_t vecmax = veclen - 1;
	const size_t freqspace = 32 * veclen;
	const size_t rarespace = 1;

	const int32_t *stopFreq = freq + lenFreq - freqspace;
	const int32_t *stopRare = rare + lenRare - rarespace;
	if (freq > stopFreq) {
		//A = freq?
		count = scalar(freq, lenFreq, rare, lenRare, out);
		return count;
	}
	for (; rare < stopRare; ++rare) {
		const int32_t matchRare = *rare; //nextRare;
		const vec Match = _mm_set1_epi32(matchRare);

		if (freq[veclen * 31 + vecmax] < matchRare) { // if no match possible
			int32_t offset = 1;
			if (freq + veclen  * 32 > stopFreq) {
				freq += veclen * 32;
				goto FINISH_SCALAR;
			}
			while (freq[veclen * offset * 32 + veclen * 31 + vecmax]
					< matchRare) { // if no match possible
				if (freq + veclen * (2 * offset) * 32 <= stopFreq) {
					offset *= 2;
				} else if (freq + veclen * (offset + 1) * 32 <= stopFreq) {
					offset = static_cast<int32_t>((stopFreq - freq) / (veclen * 32));
					//offset += 1;
					if (freq[veclen * offset * 32 + veclen * 31 + vecmax]
							< matchRare) {
						freq += veclen * offset * 32;
						goto FINISH_SCALAR;
					} else {
						break;
					}
				} else {
					freq += veclen * offset * 32;
					goto FINISH_SCALAR;
				}
			}
			int32_t lower = offset / 2;
			while (lower + 1 != offset) {
				const int32_t mid = (lower + offset) / 2;
				if (freq[veclen * mid * 32 + veclen * 31 + vecmax]
						< matchRare)
					lower = mid;
				else
					offset = mid;
			}
			freq += veclen * offset * 32;
		}
		vec Q0, Q1, Q2, Q3;
		if (freq[veclen * 15 + vecmax] >= matchRare) {
			if (freq[veclen * 7 + vecmax] < matchRare) {
				//there are 16 SSE registers in AVX architecture, we use 12 here + 1 + 1 + 1
				//if the compiler is good this code should just use 15 registers (could spill though)
				//these are all vectors we are comparing can get index from that (+8 = +8 vectors)
				//Match, r0-r7, Q0-Q3, lr Match (8+4+2)
				const size_t offset = 8;
				__m128i lr = _mm_loadu_si128(reinterpret_cast<const vec *>(freq) + offset);
				const __m128i r0 = _mm_cmpeq_epi32(lr, Match);
				lr = _mm_loadu_si128(reinterpret_cast<const vec *>(freq) + offset + 1);
				const __m128i r1 = _mm_cmpeq_epi32(lr, Match);
				Q0 = _mm_or_si128(r0,r1);

				lr = _mm_loadu_si128(reinterpret_cast<const vec *>(freq) + offset + 2);
				const __m128i r2 = _mm_cmpeq_epi32(lr, Match);
				lr = _mm_loadu_si128(reinterpret_cast<const vec *>(freq) + offset + 3);
				const __m128i r3 = _mm_cmpeq_epi32(lr, Match);
				Q1 = _mm_or_si128(r2,r3);

				lr = _mm_loadu_si128(reinterpret_cast<const vec *>(freq) + offset + 4);
				const __m128i r4 = _mm_cmpeq_epi32(lr, Match);
				lr = _mm_loadu_si128(reinterpret_cast<const vec *>(freq) + offset + 5);
				const __m128i r5 = _mm_cmpeq_epi32(lr, Match);
				Q2 = _mm_or_si128(r4,r5);

				lr = _mm_loadu_si128(reinterpret_cast<const vec *>(freq) + offset + 6);
				const __m128i r6 = _mm_cmpeq_epi32(lr, Match);
				lr = _mm_loadu_si128(reinterpret_cast<const vec *>(freq) + offset + 7);
				const __m128i r7 = _mm_cmpeq_epi32(lr, Match);
				Q3 = _mm_or_si128(r6,r7);

				lr = _mm_or_si128(Q0, Q1);
				const vec F0 = _mm_or_si128(lr, _mm_or_si128(Q2, Q3));

				if (_mm_testz_si128(F0, F0) == 0) {
					//const size_t freqOffset = offset*4 + N::check_registers(Q0,Q1,Q2,Q3,r0,r1,r2,r3,r4,r5,r6,r7);
					//const size_t hit_amount = scalar(matchRare, out, f, (rare-startRare), ((freq+freqOffset)-startFreq));
					//out = N::advanceC(out,hit_amount);
					//count += hit_amount;
					out[count++] = matchRare;
				} 
			} else {
				//there are 16 SSE registers in AVX architecture, we use 12 here + 1 + 1 + 1
				//if the compiler is good this code should just use 15 registers (could spill though)
				//these are all vectors we are comparing can get index from that (+8 = +8 vectors)
				//Match, r0-r7, Q0-Q3, lr Match (8+4+2)
				const size_t offset = 0;
				__m128i lr = _mm_loadu_si128(reinterpret_cast<const vec *>(freq) + offset);
				const __m128i r0 = _mm_cmpeq_epi32(lr, Match);
				lr = _mm_loadu_si128(reinterpret_cast<const vec *>(freq) + offset + 1);
				const __m128i r1 = _mm_cmpeq_epi32(lr, Match);
				Q0 = _mm_or_si128(r0,r1);

				lr = _mm_loadu_si128(reinterpret_cast<const vec *>(freq) + offset + 2);
				const __m128i r2 = _mm_cmpeq_epi32(lr, Match);
				lr = _mm_loadu_si128(reinterpret_cast<const vec *>(freq) + offset + 3);
				const __m128i r3 = _mm_cmpeq_epi32(lr, Match);
				Q1 = _mm_or_si128(r2,r3);

				lr = _mm_loadu_si128(reinterpret_cast<const vec *>(freq) + offset + 4);
				const __m128i r4 = _mm_cmpeq_epi32(lr, Match);
				lr = _mm_loadu_si128(reinterpret_cast<const vec *>(freq) + offset + 5);
				const __m128i r5 = _mm_cmpeq_epi32(lr, Match);
				Q2 = _mm_or_si128(r4,r5);

				lr = _mm_loadu_si128(reinterpret_cast<const vec *>(freq) + offset + 6);
				const __m128i r6 = _mm_cmpeq_epi32(lr, Match);
				lr = _mm_loadu_si128(reinterpret_cast<const vec *>(freq) + offset + 7);
				const __m128i r7 = _mm_cmpeq_epi32(lr, Match);
				Q3 = _mm_or_si128(r6,r7);

				lr = _mm_or_si128(Q0, Q1);
				const vec F0 = _mm_or_si128(lr, _mm_or_si128(Q2, Q3));

				if (_mm_testz_si128(F0, F0) == 0) {
					//const size_t freqOffset = offset*4 + N::check_registers(Q0,Q1,Q2,Q3,r0,r1,r2,r3,r4,r5,r6,r7);
					//const size_t hit_amount = N::scalar(matchRare,out,f,(rare-startRare),((freq+freqOffset)-startFreq));
					//out = N::advanceC(out,hit_amount);
					//count += hit_amount;
					out[count++] = matchRare;
				}  
			}
		} else {
			if (freq[veclen * 23 + vecmax] < matchRare) {
				//there are 16 SSE registers in AVX architecture, we use 12 hVertexID ere + 1 + 1 + 1
				//if the compiler is good this code should just use 15 registers (could spill though)
				//these are all vectors we are comparing can get index from that (+8 = +8 vectors)
				//Match, r0-r7, Q0-Q3, lr Match (8+4+2)
				const size_t offset = 24;
				__m128i lr = _mm_loadu_si128(reinterpret_cast<const vec *>(freq) + offset);
				const __m128i r0 = _mm_cmpeq_epi32(lr, Match);
				lr = _mm_loadu_si128(reinterpret_cast<const vec *>(freq) + offset + 1);
				const __m128i r1 = _mm_cmpeq_epi32(lr, Match);
				Q0 = _mm_or_si128(r0,r1);

				lr = _mm_loadu_si128(reinterpret_cast<const vec *>(freq) + offset + 2);
				const __m128i r2 = _mm_cmpeq_epi32(lr, Match);
				lr = _mm_loadu_si128(reinterpret_cast<const vec *>(freq) + offset + 3);
				const __m128i r3 = _mm_cmpeq_epi32(lr, Match);
				Q1 = _mm_or_si128(r2,r3);

				lr = _mm_loadu_si128(reinterpret_cast<const vec *>(freq) + offset + 4);
				const __m128i r4 = _mm_cmpeq_epi32(lr, Match);
				lr = _mm_loadu_si128(reinterpret_cast<const vec *>(freq) + offset + 5);
				const __m128i r5 = _mm_cmpeq_epi32(lr, Match);
				Q2 = _mm_or_si128(r4,r5);

				lr = _mm_loadu_si128(reinterpret_cast<const vec *>(freq) + offset + 6);
				const __m128i r6 = _mm_cmpeq_epi32(lr, Match);
				lr = _mm_loadu_si128(reinterpret_cast<const vec *>(freq) + offset + 7);
				const __m128i r7 = _mm_cmpeq_epi32(lr, Match);
				Q3 = _mm_or_si128(r6,r7);

				lr = _mm_or_si128(Q0, Q1);
				const vec F0 = _mm_or_si128(lr, _mm_or_si128(Q2, Q3));

				if (_mm_testz_si128(F0, F0) == 0) {
					//const size_t freqOffset = offset*4 + N::check_registers(Q0,Q1,Q2,Q3,r0,r1,r2,r3,r4,r5,r6,r7);
					//const size_t hit_amount = N::scalar(matchRare,out,f,(rare-startRare),((freq+freqOffset)-startFreq));
					//out = N::advanceC(out,hit_amount);
					//count += hit_amount;
					out[count++] = matchRare;
				} 
			} else {
				//there are 16 SSE registers in AVX architecture, we use 12 here + 1 + 1 + 1
				//if the compiler is good this code should just use 15 registers (could spill though)
				//these are all vectors we are comparing can get index from that (+8 = +8 vectors)
				//Match, r0-r7, Q0-Q3, lr Match (8+4+2)
				const size_t offset = 16;
				__m128i lr = _mm_loadu_si128(reinterpret_cast<const vec *>(freq) + offset);
				const __m128i r0 = _mm_cmpeq_epi32(lr, Match);
				lr = _mm_loadu_si128(reinterpret_cast<const vec *>(freq) + offset + 1);
				const __m128i r1 = _mm_cmpeq_epi32(lr, Match);
				Q0 = _mm_or_si128(r0,r1);

				lr = _mm_loadu_si128(reinterpret_cast<const vec *>(freq) + offset + 2);
				const __m128i r2 = _mm_cmpeq_epi32(lr, Match);
				lr = _mm_loadu_si128(reinterpret_cast<const vec *>(freq) + offset + 3);
				const __m128i r3 = _mm_cmpeq_epi32(lr, Match);
				Q1 = _mm_or_si128(r2,r3);

				lr = _mm_loadu_si128(reinterpret_cast<const vec *>(freq) + offset + 4);
				const __m128i r4 = _mm_cmpeq_epi32(lr, Match);
				lr = _mm_loadu_si128(reinterpret_cast<const vec *>(freq) + offset + 5);
				const __m128i r5 = _mm_cmpeq_epi32(lr, Match);
				Q2 = _mm_or_si128(r4,r5);

				lr = _mm_loadu_si128(reinterpret_cast<const vec *>(freq) + offset + 6);
				const __m128i r6 = _mm_cmpeq_epi32(lr, Match);
				lr = _mm_loadu_si128(reinterpret_cast<const vec *>(freq) + offset + 7);
				const __m128i r7 = _mm_cmpeq_epi32(lr, Match);
				Q3 = _mm_or_si128(r6,r7);

				lr = _mm_or_si128(Q0, Q1);
				const vec F0 = _mm_or_si128(lr, _mm_or_si128(Q2, Q3));

				if (_mm_testz_si128(F0, F0) == 0) {
					//const size_t freqOffset = offset*4 + N::check_registers(Q0,Q1,Q2,Q3,r0,r1,r2,r3,r4,r5,r6,r7);
					//const size_t hit_amount = N::scalar(matchRare,out,f,(rare-startRare),((freq+freqOffset)-startFreq));
					//out = N::advanceC(out,hit_amount);
					//count += hit_amount;
					out[count++] = matchRare;
				} 
			}

		}
	}

FINISH_SCALAR: 
	//A = rare?
	size_t final_count = count + scalar(rare, stopRare + rarespace - rare, freq, stopFreq + freqspace - freq, out + count);
#ifdef PRINT_INTERSECTION
    if (calculate_domains_context)
        CD_totalCountGalloping += final_count;
    else
        totalCountGalloping += final_count;
#endif
	return final_count;
#endif
}

int intersect_shuffling(const int32_t* A, const int32_t* B, size_t s_a, size_t s_b, int* C) 
{
#ifdef PRINT_INTERSECTION
    if (calculate_domains_context) {
        CD_totalA += s_a;
        CD_totalB += s_b;
        CD_totalShuffling++;
    }
    else {
        totalA += s_a;
        totalB += s_b;
        totalShuffling++;
    }
#endif
	size_t count = 0;
	size_t i_a = 0, i_b = 0;

#ifdef OPT_INTERSECT
	if(A[0] < B[0]){
		int first = 0, last = s_a;
		while(first != last)
		{
			int mid = (first + last) / 2;
			if(A[mid] < B[0])
                first = mid + 1;
			else
                last = mid;	
		}
		i_a = (last / 4) * 4;
	}
	else if(A[0] > B[0]){
		int first = 0, last = s_b;
		while(first != last)
		{
			int mid = (first + last) / 2;
			if(B[mid] < A[0])
                first = mid + 1;
			else
                last = mid;	
		}
		i_b = (last / 4) * 4;
	}
#endif
    
#ifdef SIMD_SHUFFLING
    //assert(false);
	size_t st_a = (s_a / 4) * 4;
    size_t st_b = (s_b / 4) * 4;

	while(i_a < st_a && i_b < st_b){
		__m128i v_a = _mm_loadu_si128((__m128i*)&A[i_a]);
		__m128i v_b = _mm_loadu_si128((__m128i*)&B[i_b]);

		VertexID a_max = _mm_extract_epi32(v_a, 3);
		VertexID b_max = _mm_extract_epi32(v_b, 3);
#ifdef OPT_INTERSECT
		VertexID a_min = _mm_extract_epi32(v_a, 0);
		VertexID b_min = _mm_extract_epi32(v_b, 0);

		if(a_max < b_min){
			i_a += 4;
			continue;
		}
		if(b_max < a_min){
			i_b += 4;
			continue;
		}
#endif
		i_a += (a_max <= b_max) * 4;
		i_b += (a_max >= b_max) * 4;

		int32_t cyclic_shift = _MM_SHUFFLE(0,3,2,1);
		__m128i cmp_mask1 = _mm_cmpeq_epi32(v_a, v_b);    // pairwise comparison
		v_b = _mm_shuffle_epi32(v_b, cyclic_shift);       // shuffling
		__m128i cmp_mask2 = _mm_cmpeq_epi32(v_a, v_b);    // again...
		v_b = _mm_shuffle_epi32(v_b, cyclic_shift);
		__m128i cmp_mask3 = _mm_cmpeq_epi32(v_a, v_b);    // and again...
		v_b = _mm_shuffle_epi32(v_b, cyclic_shift);
		__m128i cmp_mask4 = _mm_cmpeq_epi32(v_a, v_b);    // and again.
		__m128i cmp_mask = _mm_or_si128(
				_mm_or_si128(cmp_mask1, cmp_mask2),
				_mm_or_si128(cmp_mask3, cmp_mask4)
				); // OR-ing of comparison masks
		// convert the 128-bit mask to the 4-bit mask
		int32_t mask = _mm_movemask_ps((__m128)cmp_mask);

		// copy out common elements
		__m128i p = _mm_shuffle_epi8(v_a, shuffle_mask[mask]);
		_mm_storeu_si128((__m128i*)&C[count], p);
		count += _mm_popcnt_u32(mask); // a number of elements is a weight of the mask
	}
#endif
    
#ifdef SIMPLE_SCALAR
    //this is faster than scalar function..! why?
    while(i_a < s_a && i_b < s_b)
    {
        if(A[i_a] < B[i_b])
            i_a++;
        else if(B[i_b] < A[i_a])
            i_b++;
        else {
            C[count++] = A[i_a];
            i_a++;
            i_b++;
        }
    }
#else
    count += scalar(A+i_a, s_a-i_a, B+i_b, s_b-i_b, C+count);
#endif

#ifdef PRINT_INTERSECTION
    if (calculate_domains_context)
        CD_totalCountShuffling += count;
    else
        totalCountShuffling += count;
#endif
	return count;
}

#endif

/***************************** API *****************************/

#ifdef COMPRESS
//binary
int intersect(range& RA, range& RB, int* CB, int* CS) //, int& NB, int& NV)
{
	const int* A = RA.bases.first;
	const int* B = RB.bases.first;

	size_t s_a = RA.bases.second - A; 
	size_t s_b = RB.bases.second - B; 

	if(s_a == 0 || s_b == 0)
    {
        //NB = NV = 0;
        return 0;
    }
	else if(A[s_a-1] < B[0] || B[s_b-1] < A[0]) //20 ms
    {
        //NB = NV = 0;
        return 0;
    }
    else
    {
        if(s_a < s_b){
            if(32*s_a > s_b)
            {
                return intersect_qfilter_bsr_b4_v2(
                        A, RA.states.first, s_a,
                        B, RB.states.first, s_b,
                        CB, CS); //, NB, NV);
            }
            else
            {
                return intersect_simdgalloping_bsr(
                        A, RA.states.first, s_a,
                        B, RB.states.first, s_b,
                        CB, CS); //, NB, NV);
            }
        }
        else{
            if(32*s_b > s_a)
            {
                return intersect_qfilter_bsr_b4_v2(
                        B, RB.states.first, s_b,
                        A, RA.states.first, s_a,
                        CB, CS); //, NB, NV);
            }
            else
            {
                return intersect_simdgalloping_bsr(
                        B, RB.states.first, s_b,
                        A, RA.states.first, s_a,
                        CB, CS); //, NB, NV);
            }
        }
    }
}

//assume that R has >= 2 ranges
int intersect(std::vector<range>& R, int* CB, int* CS, int* TB, int* TS) //, int& NB, int& NV)
{
    int count = 0;
    int* next_B;
    int* next_S;

    if((R.size() % 2) == 1)
    {
		//switch temp and C
        next_B = CB;
        next_S = CS;
        CB = TB;
        CS = TS;
        TB = next_B;
        TS = next_S;
    }
	//intersect first two
    int s_a = R[0].bases.second - R[0].bases.first;
    int s_b = R[1].bases.second - R[1].bases.first;
#ifdef VERBOSE
    std::cout << "A # base = " << s_a << ", B # base = " << s_b << std::endl;
#endif
	count = intersect(R[0].bases.first, R[0].states.first, s_a,
			          R[1].bases.first, R[1].states.first, s_b, CB, CS);

	int i = 2;
	//then intersect one at a time
	while (i < R.size()) {
		//switch temp and C
        next_B = CB;
        next_S = CS;
        CB = TB;
        CS = TS;
        TB = next_B;
        TS = next_S;
#ifdef VERBOSE
        std::cout << "C # base = " << count << ", next # base = " << (iter->bases.second - iter->bases.first) << std::endl;
#endif
		count = intersect(TB, TS, count, 
				R[i].bases.first, R[i].states.first, R[i].bases.second - R[i].bases.first,
				CB, CS);
		i++;
	}
#ifdef VERBOSE
    std::cout << "final C # base = " << count << std::endl;
#endif
    return count;
}

#else 

//bottleneck!
int intersect(const range& RA, const range& RB, int* C) 
{
	const int32_t * A = RA.begin;
	const int32_t * B = RB.begin;

	size_t s_a = RA.end - A; 
	size_t s_b = RB.end - B; 

	if(s_a == 0 || s_b == 0) {
#ifdef PRINT_INTERSECTION
        if (calculate_domains_context)
            CD_totalEdgeCase1++;
        else
            totalEdgeCase1++;
#endif
		return 0;
    }
	if(A[s_a-1] < B[0] || B[s_b-1] < A[0]) {
#ifdef PRINT_INTERSECTION
        if (calculate_domains_context)
            CD_totalEdgeCase2++;
        else
            totalEdgeCase2++;
#endif
		return 0;
    }

	if(s_a < s_b){
		if(32*s_a < s_b)
			return intersect_galloping(A, B, s_a, s_b, C);
		else
			return intersect_shuffling(A, B, s_a, s_b, C);
	}
	else{
		if(32*s_b < s_a)
			return intersect_galloping(B, A, s_b, s_a, C);
		else
			return intersect_shuffling(B, A, s_b, s_a, C);
	}
}

int intersect(const meta_range& RA, const meta_range& RB, int* C) 
{
	const int32_t * A = RA.begin;
	const int32_t * B = RB.begin;

	size_t s_a = RA.end - A; 
	size_t s_b = RB.end - B; 

	if(s_a == 0 || s_b == 0) {
		return 0;
    }
	if(A[s_a-1] < B[0] || B[s_b-1] < A[0]) {
		return 0;
    }

	if(s_a < s_b){
		if(32*s_a < s_b)
			return intersect_galloping(A, B, s_a, s_b, C);
		else
			return intersect_shuffling(A, B, s_a, s_b, C);
	}
	else{
		if(32*s_b < s_a)
			return intersect_galloping(B, A, s_b, s_a, C);
		else
			return intersect_shuffling(B, A, s_b, s_a, C);
	}
}

//assume S < RB size
int intersect(int *V, size_t S, const range& RB, int* C)
{
	const int32_t * B = RB.begin;

	size_t s_b = RB.end - B; 

	if(s_b == 0) {
#ifdef PRINT_INTERSECTION
        if (calculate_domains_context)
            CD_totalEdgeCase1++;
        else
            totalEdgeCase1++;
#endif
		return 0;
    }
	if(V[S-1] < B[0] || B[s_b-1] < V[0]) {
#ifdef PRINT_INTERSECTION
        if (calculate_domains_context)
            CD_totalEdgeCase2++;
        else
            totalEdgeCase2++;
#endif
		return 0;
    }

	if(32*S < s_b)
		return intersect_galloping(V, B, S, s_b, C);
	else
		return intersect_shuffling(V, B, S, s_b, C);
}

int intersect(int *V, size_t S, const meta_range& RB, int* C)
{
	const int32_t * B = RB.begin;

	size_t s_b = RB.end - B; 

	if(s_b == 0) {
		return 0;
    }
	if(V[S-1] < B[0] || B[s_b-1] < V[0]) {
		return 0;
    }

	if(32*S < s_b)
		return intersect_galloping(V, B, S, s_b, C);
	else
		return intersect_shuffling(V, B, S, s_b, C);
}

#ifdef OPT
int intersect(std::vector<range>& R, int* C, int* temp)
#else
int intersect(std::set<range>& R, int* C, int* temp)
#endif
{
	int count = 0;
	int *next;
    if((R.size() % 2) == 1)
    {
		//switch temp and C
        next = C;
        C = temp;
        temp = next;
    }
	//intersect first two
#ifdef OPT
    count = intersect(R[0], R[1], C);
	int i = 2;
#else
	std::set<range>::iterator t = R.begin();
	std::set<range>::iterator iter = ++t; 
	std::set<range>::iterator first = R.begin();
	count = intersect(*first, *iter, C);
	iter++;
#endif

    if (count == 0)
        return 0;
		
	//then intersect one at a time
#ifdef OPT
	while (i < R.size()) {
#else
	while (iter != R.end()) {
#endif
		//switch temp and C
		next = C;
		C = temp;
		temp = next;
#ifdef OPT
		count = intersect(temp, count, R[i], C);
		i++;
#else
		count = intersect(temp, count, *iter, C);
		iter++;
#endif
        if (count == 0)
            return 0;
	}

	return count;
}

int intersect(std::set<meta_range>& R, int* C, int* temp, AdjElem& failed_edge)
{
	int count = 0;
	int *next;
    if((R.size() % 2) == 1)
    {
		//switch temp and C
        next = C;
        C = temp;
        temp = next;
    }
	//intersect first two

	std::set<meta_range>::iterator t = R.begin();
	std::set<meta_range>::iterator iter = ++t; 
	std::set<meta_range>::iterator first = R.begin();
	count = intersect(*first, *iter, C);
	if (count == 0) {
		failed_edge = iter->edge;
		return count;
	}
	iter++;
		
	//then intersect one at a time
	while (iter != R.end()) {
		//switch temp and C
		next = C;
		C = temp;
		temp = next;

		count = intersect(temp, count, *iter, C);
		if (count == 0) {
			failed_edge = iter->edge;
			return count;
		}
		iter++;
	}

	return count;
}

int do_union(const range& RA, const range& RB, int* C) 
{
	const int32_t * A = RA.begin;
	const int32_t * B = RB.begin;

	size_t s_a = RA.end - A; 
	size_t s_b = RB.end - B; 

	if(s_a == 0 || s_b == 0)
		return 0;

    int i = 0, j = 0, k = 0;
    while (i < s_a && j < s_b) {
        if (A[i] < B[j])
            C[k++] = A[i++];
        else
            C[k++] = B[j++];
    }
    while (i < s_a)
        C[k++] = A[i++];
    while (j < s_b)
        C[k++] = B[j++];

    return k;
}

int do_union(int *V, size_t S, const range& RB, int* C)
{
	const int32_t * B = RB.begin;

	size_t s_b = RB.end - B; 

	if(S == 0 || s_b == 0)
		return 0;

    int i = 0, j = 0, k = 0;
    while (i < S && j < s_b) {
        if (V[i] < B[j])
            C[k++] = V[i++];
        else
            C[k++] = B[j++];
    }
    while (i < S)
        C[k++] = V[i++];
    while (j < s_b)
        C[k++] = B[j++];

    return k;
}

int do_union(std::vector<range>& ranges, int* C)
{
	//vid, range index
    typedef pair<int, int> elem;
    priority_queue<elem, vector<elem>, greater<elem>> q;
    for (int j = 0; j < ranges.size(); j++) {
		assert(ranges[j].end != ranges[j].begin);
        q.emplace(ranges[j].begin[0], j);
        ranges[j].begin++;
    }
    int k = 0;
    int cur_id = -1;
    while (!q.empty()) {
        pair<int, int> t = q.top();
        q.pop();
        if (t.first > cur_id) {
            cur_id = t.first;
            C[k++] = cur_id;
        }
        if (ranges[t.second].end != ranges[t.second].begin) {
            q.emplace(*ranges[t.second].begin, t.second); 
            ranges[t.second].begin++;
        }
    }

	return k;
}

int do_union(std::set<range>& ranges, int* C)
{
	std::vector<range> vec(ranges.begin(), ranges.end());
	return do_union(vec, C);
}

#endif
