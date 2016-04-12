#line 2 "luxrays_types.cl"

/***************************************************************************
 * Copyright 1998-2015 by authors (see AUTHORS.txt)                        *
 *                                                                         *
 *   This file is part of LuxRender.                                       *
 *                                                                         *
 * Licensed under the Apache License, Version 2.0 (the "License");         *
 * you may not use this file except in compliance with the License.        *
 * You may obtain a copy of the License at                                 *
 *                                                                         *
 *     http://www.apache.org/licenses/LICENSE-2.0                          *
 *                                                                         *
 * Unless required by applicable law or agreed to in writing, software     *
 * distributed under the License is distributed on an "AS IS" BASIS,       *
 * WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.*
 * See the License for the specific language governing permissions and     *
 * limitations under the License.                                          *
 ***************************************************************************/

#define NULL_INDEX (0xffffffffu)

#if defined(LUXRAYS_OPENCL_KERNEL)

#define NULL 0

#if defined(__APPLE_CL__)
float3 __OVERLOAD__ mix(float3 a, float3 b, float t)
{
	return a + ( b - a ) * t;
}
#endif

#if defined(__APPLE_FIX__)

float2 VLOAD2F(const __global float *p) {
	return (float2)(p[0], p[1]);
}

void VSTORE2F(const float2 v, __global float *p) {
	p[0] = v.x;
	p[1] = v.y;
}

float3 VLOAD3F(const __global float *p) {
	return (float3)(p[0], p[1], p[2]);
}

float3 VLOAD3F_Private(const float *p) {
	return (float3)(p[0], p[1], p[2]);
}

void VSTORE3F(const float3 v, __global float *p) {
	p[0] = v.x;
	p[1] = v.y;
	p[2] = v.z;
}

float4 VLOAD4F(const __global float *p) {
	return (float4)(p[0], p[1], p[2], p[3]);
}

float4 VLOAD4F_Private(const float *p) {
	return (float4)(p[0], p[1], p[2], p[3]);
}

void VSTORE4F(const float4 v, __global float *p) {
	p[0] = v.x;
	p[1] = v.y;
	p[2] = v.z;
	p[3] = v.w;
}

#else

float2 VLOAD2F(const __global float *p) {
	return vload2(0, p);
}

void VSTORE2F(const float2 v, __global float *p) {
	vstore2(v, 0, p);
}

float3 VLOAD3F(const __global float *p) {
	return vload3(0, p);
}

float3 VLOAD3F_Private(const float *p) {
	return vload3(0, p);
}

void VSTORE3F(const float3 v, __global float *p) {
	vstore3(v, 0, p);
}

float4 VLOAD4F(const __global float *p) {
	return vload4(0, p);
}

// Input address must be aligned to 16B
// This performs better than vload4()
float4 VLOAD4F_Align(const __global float *p) {
	return *((const __global float4 *)p);
}

float4 VLOAD4F_Private(const float *p) {
	return vload4(0, p);
}

void VSTORE4F(const float4 v, __global float *p) {
	vstore4(v, 0, p);
}

#endif

void VADD3F(__global float *p, const float3 v) {
	VSTORE3F(VLOAD3F(p) + v, p);
}

#endif
#line 2 "epsilon_types.cl"

/***************************************************************************
 * Copyright 1998-2015 by authors (see AUTHORS.txt)                        *
 *                                                                         *
 *   This file is part of LuxRender.                                       *
 *                                                                         *
 * Licensed under the Apache License, Version 2.0 (the "License");         *
 * you may not use this file except in compliance with the License.        *
 * You may obtain a copy of the License at                                 *
 *                                                                         *
 *     http://www.apache.org/licenses/LICENSE-2.0                          *
 *                                                                         *
 * Unless required by applicable law or agreed to in writing, software     *
 * distributed under the License is distributed on an "AS IS" BASIS,       *
 * WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.*
 * See the License for the specific language governing permissions and     *
 * limitations under the License.                                          *
 ***************************************************************************/

// NOTE: DEFAULT_EPSILON_MIN is very small. A plane passing exactly for the
// origin will suffer of self shadow problems because the Ray class will use
// MachineEpsilon(ray.o) as epsilon for the ray.mint. However it is pretty much
// the only case where there is a problem so better to not change anything.
// As workaround, moving the plane away from the origin is enough.
#define DEFAULT_EPSILON_MIN 1e-9f
#define DEFAULT_EPSILON_MAX 1e-1f
#define DEFAULT_EPSILON_STATIC 1e-5f

// An epsilon that can be used as threshold for cos(theta). For instance:
// if (Dot(N, LightDir) < DEFAULT_COS_EPSILON_STATIC) return Spectrum();
#define DEFAULT_COS_EPSILON_STATIC 1e-4f

// This is about 1e-5f for values near 1.f
#define DEFAULT_EPSILON_DISTANCE_FROM_VALUE 0x80u
#line 2 "epsilon_funcs.cl"

/***************************************************************************
 * Copyright 1998-2015 by authors (see AUTHORS.txt)                        *
 *                                                                         *
 *   This file is part of LuxRender.                                       *
 *                                                                         *
 * Licensed under the Apache License, Version 2.0 (the "License");         *
 * you may not use this file except in compliance with the License.        *
 * You may obtain a copy of the License at                                 *
 *                                                                         *
 *     http://www.apache.org/licenses/LICENSE-2.0                          *
 *                                                                         *
 * Unless required by applicable law or agreed to in writing, software     *
 * distributed under the License is distributed on an "AS IS" BASIS,       *
 * WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.*
 * See the License for the specific language governing permissions and     *
 * limitations under the License.                                          *
 ***************************************************************************/

float MachineEpsilon_FloatAdvance(const float value) {
	return as_float(as_uint(value) + DEFAULT_EPSILON_DISTANCE_FROM_VALUE);
}

float MachineEpsilon_E(const float value) {
	const float epsilon = fabs(MachineEpsilon_FloatAdvance(value) - value);

	return clamp(epsilon, PARAM_RAY_EPSILON_MIN, PARAM_RAY_EPSILON_MAX);
}

float MachineEpsilon_E_Float3(const float3 v) {
	return fmax(MachineEpsilon_E(v.x), fmax(MachineEpsilon_E(v.y), MachineEpsilon_E(v.z)));
}
#line 2 "point_types.cl"

/***************************************************************************
 * Copyright 1998-2015 by authors (see AUTHORS.txt)                        *
 *                                                                         *
 *   This file is part of LuxRender.                                       *
 *                                                                         *
 * Licensed under the Apache License, Version 2.0 (the "License");         *
 * you may not use this file except in compliance with the License.        *
 * You may obtain a copy of the License at                                 *
 *                                                                         *
 *     http://www.apache.org/licenses/LICENSE-2.0                          *
 *                                                                         *
 * Unless required by applicable law or agreed to in writing, software     *
 * distributed under the License is distributed on an "AS IS" BASIS,       *
 * WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.*
 * See the License for the specific language governing permissions and     *
 * limitations under the License.                                          *
 ***************************************************************************/

typedef struct {
	float x, y, z;
} Point;
#line 2 "vector_types.cl"

/***************************************************************************
 * Copyright 1998-2015 by authors (see AUTHORS.txt)                        *
 *                                                                         *
 *   This file is part of LuxRender.                                       *
 *                                                                         *
 * Licensed under the Apache License, Version 2.0 (the "License");         *
 * you may not use this file except in compliance with the License.        *
 * You may obtain a copy of the License at                                 *
 *                                                                         *
 *     http://www.apache.org/licenses/LICENSE-2.0                          *
 *                                                                         *
 * Unless required by applicable law or agreed to in writing, software     *
 * distributed under the License is distributed on an "AS IS" BASIS,       *
 * WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.*
 * See the License for the specific language governing permissions and     *
 * limitations under the License.                                          *
 ***************************************************************************/

#define ASSIGN_VECTOR(a, b) { (a).x = (b).x; (a).y = (b).y; (a).z = (b).z; }

typedef struct {
	float x, y, z;
} Vector;
#line 2 "matrix4x4_types.cl"

/***************************************************************************
 * Copyright 1998-2015 by authors (see AUTHORS.txt)                        *
 *                                                                         *
 *   This file is part of LuxRender.                                       *
 *                                                                         *
 * Licensed under the Apache License, Version 2.0 (the "License");         *
 * you may not use this file except in compliance with the License.        *
 * You may obtain a copy of the License at                                 *
 *                                                                         *
 *     http://www.apache.org/licenses/LICENSE-2.0                          *
 *                                                                         *
 * Unless required by applicable law or agreed to in writing, software     *
 * distributed under the License is distributed on an "AS IS" BASIS,       *
 * WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.*
 * See the License for the specific language governing permissions and     *
 * limitations under the License.                                          *
 ***************************************************************************/

typedef struct {
	float m[4][4];
} Matrix4x4;
#line 2 "matrix4x4_funcs.cl"

/***************************************************************************
 * Copyright 1998-2015 by authors (see AUTHORS.txt)                        *
 *                                                                         *
 *   This file is part of LuxRender.                                       *
 *                                                                         *
 * Licensed under the Apache License, Version 2.0 (the "License");         *
 * you may not use this file except in compliance with the License.        *
 * You may obtain a copy of the License at                                 *
 *                                                                         *
 *     http://www.apache.org/licenses/LICENSE-2.0                          *
 *                                                                         *
 * Unless required by applicable law or agreed to in writing, software     *
 * distributed under the License is distributed on an "AS IS" BASIS,       *
 * WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.*
 * See the License for the specific language governing permissions and     *
 * limitations under the License.                                          *
 ***************************************************************************/

float3 Matrix4x4_ApplyPoint(__global const Matrix4x4* restrict m, const float3 point) {
	const float4 point4 = (float4)(point.x, point.y, point.z, 1.f);

	const float4 row3 = VLOAD4F(&m->m[3][0]);
	const float iw = 1.f / dot(row3, point4);

	const float4 row0 = VLOAD4F(&m->m[0][0]);
	const float4 row1 = VLOAD4F(&m->m[1][0]);
	const float4 row2 = VLOAD4F(&m->m[2][0]);
	return (float3)(
			iw * dot(row0, point4),
			iw * dot(row1, point4),
			iw * dot(row2, point4)
			);
}

float3 Matrix4x4_ApplyPoint_Align(__global const Matrix4x4* restrict m, const float3 point) {
	const float4 point4 = (float4)(point.x, point.y, point.z, 1.f);

	const float4 row3 = VLOAD4F_Align(&m->m[3][0]);
	const float iw = 1.f / dot(row3, point4);

	const float4 row0 = VLOAD4F_Align(&m->m[0][0]);
	const float4 row1 = VLOAD4F_Align(&m->m[1][0]);
	const float4 row2 = VLOAD4F_Align(&m->m[2][0]);
	return (float3)(
			iw * dot(row0, point4),
			iw * dot(row1, point4),
			iw * dot(row2, point4)
			);
}

float3 Matrix4x4_ApplyPoint_Private(Matrix4x4 *m, const float3 point) {
	const float4 point4 = (float4)(point.x, point.y, point.z, 1.f);

	const float4 row3 = VLOAD4F_Private(&m->m[3][0]);
	const float iw = 1.f / dot(row3, point4);

	const float4 row0 = VLOAD4F_Private(&m->m[0][0]);
	const float4 row1 = VLOAD4F_Private(&m->m[1][0]);
	const float4 row2 = VLOAD4F_Private(&m->m[2][0]);
	return (float3)(
			iw * dot(row0, point4),
			iw * dot(row1, point4),
			iw * dot(row2, point4)
			);
}

float3 Matrix4x4_ApplyVector(__global const Matrix4x4* restrict m, const float3 vector) {
	const float3 row0 = VLOAD3F(&m->m[0][0]);
	const float3 row1 = VLOAD3F(&m->m[1][0]);
	const float3 row2 = VLOAD3F(&m->m[2][0]);
	return (float3)(
			dot(row0, vector),
			dot(row1, vector),
			dot(row2, vector)
			);
}

float3 Matrix4x4_ApplyVector_Private(Matrix4x4 *m, const float3 vector) {
	const float3 row0 = VLOAD3F_Private(&m->m[0][0]);
	const float3 row1 = VLOAD3F_Private(&m->m[1][0]);
	const float3 row2 = VLOAD3F_Private(&m->m[2][0]);
	return (float3)(
			dot(row0, vector),
			dot(row1, vector),
			dot(row2, vector)
			);
}

float3 Matrix4x4_ApplyNormal(__global const Matrix4x4* restrict m, const float3 normal) {
	const float3 row0 = (float3)(m->m[0][0], m->m[1][0], m->m[2][0]);
	const float3 row1 = (float3)(m->m[0][1], m->m[1][1], m->m[2][1]);
	const float3 row2 = (float3)(m->m[0][2], m->m[1][2], m->m[2][2]);
	return (float3)(
			dot(row0, normal),
			dot(row1, normal),
			dot(row2, normal)
			);
}

void Matrix4x4_Identity(Matrix4x4 *m) {
	for (int j = 0; j < 4; ++j)
		for (int i = 0; i < 4; ++i)
			m->m[i][j] = (i == j) ? 1.f : 0.f;
}

void Matrix4x4_Invert(Matrix4x4 *m) {
	int indxc[4], indxr[4];
	int ipiv[4] = {0, 0, 0, 0};

	for (int i = 0; i < 4; ++i) {
		int irow = -1, icol = -1;
		float big = 0.;
		// Choose pivot
		for (int j = 0; j < 4; ++j) {
			if (ipiv[j] != 1) {
				for (int k = 0; k < 4; ++k) {
					if (ipiv[k] == 0) {
						if (fabs(m->m[j][k]) >= big) {
							big = fabs(m->m[j][k]);
							irow = j;
							icol = k;
						}
					} else if (ipiv[k] > 1) {
						//throw std::runtime_error("Singular matrix in MatrixInvert: " + ToString(*this));
						// I can not do very much here
						Matrix4x4_Identity(m);
						return;
					}
				}
			}
		}
		++ipiv[icol];
		// Swap rows _irow_ and _icol_ for pivot
		if (irow != icol) {
			for (int k = 0; k < 4; ++k) {
				const float tmp = m->m[irow][k];
				m->m[irow][k] = m->m[icol][k];
				m->m[icol][k] = tmp;
			}
		}
		indxr[i] = irow;
		indxc[i] = icol;
		if (m->m[icol][icol] == 0.f) {
			//throw std::runtime_error("Singular matrix in MatrixInvert: " + ToString(*this));
			// I can not do very much here
			Matrix4x4_Identity(m);
			return;
		}
		// Set $m[icol][icol]$ to one by scaling row _icol_ appropriately
		float pivinv = 1.f / m->m[icol][icol];
		m->m[icol][icol] = 1.f;
		for (int j = 0; j < 4; ++j)
			m->m[icol][j] *= pivinv;
		// Subtract this row from others to zero out their columns
		for (int j = 0; j < 4; ++j) {
			if (j != icol) {
				float save = m->m[j][icol];
				m->m[j][icol] = 0;
				for (int k = 0; k < 4; ++k)
					m->m[j][k] -= m->m[icol][k] * save;
			}
		}
	}
	// Swap columns to reflect permutation
	for (int j = 3; j >= 0; --j) {
		if (indxr[j] != indxc[j]) {
			for (int k = 0; k < 4; ++k) {
				const float tmp = m->m[k][indxr[j]];
				m->m[k][indxr[j]] = m->m[k][indxc[j]];
				m->m[k][indxc[j]] = tmp;
			}
		}
	}
}
#line 2 "quaternion_types.cl"

/***************************************************************************
 * Copyright 1998-2015 by authors (see AUTHORS.txt)                        *
 *                                                                         *
 *   This file is part of LuxRender.                                       *
 *                                                                         *
 * Licensed under the Apache License, Version 2.0 (the "License");         *
 * you may not use this file except in compliance with the License.        *
 * You may obtain a copy of the License at                                 *
 *                                                                         *
 *     http://www.apache.org/licenses/LICENSE-2.0                          *
 *                                                                         *
 * Unless required by applicable law or agreed to in writing, software     *
 * distributed under the License is distributed on an "AS IS" BASIS,       *
 * WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.*
 * See the License for the specific language governing permissions and     *
 * limitations under the License.                                          *
 ***************************************************************************/

typedef struct {
	float w;
	Vector v;
} Quaternion;
#line 2 "quaternion_funcs.cl"

/***************************************************************************
 * Copyright 1998-2015 by authors (see AUTHORS.txt)                        *
 *                                                                         *
 *   This file is part of LuxRender.                                       *
 *                                                                         *
 * Licensed under the Apache License, Version 2.0 (the "License");         *
 * you may not use this file except in compliance with the License.        *
 * You may obtain a copy of the License at                                 *
 *                                                                         *
 *     http://www.apache.org/licenses/LICENSE-2.0                          *
 *                                                                         *
 * Unless required by applicable law or agreed to in writing, software     *
 * distributed under the License is distributed on an "AS IS" BASIS,       *
 * WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.*
 * See the License for the specific language governing permissions and     *
 * limitations under the License.                                          *
 ***************************************************************************/

// Get the rotation matrix from quaternion
void Quaternion_ToMatrix(const float4 q, Matrix4x4 *m) {
	const float xx = q.s1 * q.s1;
	const float yy = q.s2 * q.s2;
	const float zz = q.s3 * q.s3;
	const float xy = q.s1 * q.s2;
	const float xz = q.s1 * q.s3;
	const float yz = q.s2 * q.s3;
	const float xw = q.s1 * q.s0;
	const float yw = q.s2 * q.s0;
	const float zw = q.s3 * q.s0;

	m->m[0][0] = 1.f - 2.f * (yy + zz);
	m->m[1][0] = 2.f * (xy - zw);
	m->m[2][0] = 2.f * (xz + yw);
	m->m[0][1] = 2.f * (xy + zw);
	m->m[1][1] = 1.f - 2.f * (xx + zz);
	m->m[2][1] = 2.f * (yz - xw);
	m->m[0][2] = 2.f * (xz - yw);
	m->m[1][2] = 2.f * (yz + xw);
	m->m[2][2] = 1.f - 2.f * (xx + yy);

	// Complete matrix
	m->m[0][3] = m->m[1][3] = m->m[2][3] = 0.f;
	m->m[3][0] = m->m[3][1] = m->m[3][2] = 0.f;
	m->m[3][3] = 1.f;
}

float4 Quaternion_Slerp(float t, const float4 q1, const float4 q2) {

	float cos_phi = dot(q1, q2);
	const float sign = (cos_phi > 0.f) ? 1.f : -1.f;
	
	cos_phi *= sign;

	float f1, f2;
	if (1.f - cos_phi > 1e-6f) {	
		float phi = acos(cos_phi);
		float sin_phi = sin(phi);	
		f1 = sin((1.f - t) * phi) / sin_phi;
		f2 = sin(t * phi) / sin_phi;
	} else {
		// start and end are very close
		// perform linear interpolation
		f1 = 1.f - t;
		f2 = t;
	}

	return f1 * q1 + (sign * f2) * q2;
}
#line 2 "ray_types.cl"

/***************************************************************************
 * Copyright 1998-2015 by authors (see AUTHORS.txt)                        *
 *                                                                         *
 *   This file is part of LuxRender.                                       *
 *                                                                         *
 * Licensed under the Apache License, Version 2.0 (the "License");         *
 * you may not use this file except in compliance with the License.        *
 * You may obtain a copy of the License at                                 *
 *                                                                         *
 *     http://www.apache.org/licenses/LICENSE-2.0                          *
 *                                                                         *
 * Unless required by applicable law or agreed to in writing, software     *
 * distributed under the License is distributed on an "AS IS" BASIS,       *
 * WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.*
 * See the License for the specific language governing permissions and     *
 * limitations under the License.                                          *
 ***************************************************************************/

typedef struct {
	Point o;
	Vector d;
	float mint, maxt, time;
	float pad[3];
} Ray;

typedef struct {
	float t, b1, b2;
	unsigned int meshIndex, triangleIndex;
} RayHit;
#line 2 "ray_funcs.cl"

/***************************************************************************
 * Copyright 1998-2015 by authors (see AUTHORS.txt)                        *
 *                                                                         *
 *   This file is part of LuxRender.                                       *
 *                                                                         *
 * Licensed under the Apache License, Version 2.0 (the "License");         *
 * you may not use this file except in compliance with the License.        *
 * You may obtain a copy of the License at                                 *
 *                                                                         *
 *     http://www.apache.org/licenses/LICENSE-2.0                          *
 *                                                                         *
 * Unless required by applicable law or agreed to in writing, software     *
 * distributed under the License is distributed on an "AS IS" BASIS,       *
 * WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.*
 * See the License for the specific language governing permissions and     *
 * limitations under the License.                                          *
 ***************************************************************************/
void Ray_Init4_Private(Ray *ray, const float3 orig, const float3 dir,
		const float mint, const float maxt, const float time) {
	ray->o.x = orig.x;
	ray->o.y = orig.y;
	ray->o.z = orig.z;

	ray->d.x = dir.x;
	ray->d.y = dir.y;
	ray->d.z = dir.z;

	ray->mint = mint;
	ray->maxt = maxt;

	ray->time = time;
}

void Ray_Init3_Private(Ray *ray, const float3 orig, const float3 dir,
		const float maxt, const float time) {
	ray->o.x = orig.x;
	ray->o.y = orig.y;
	ray->o.z = orig.z;

	ray->d.x = dir.x;
	ray->d.y = dir.y;
	ray->d.z = dir.z;

	ray->mint = MachineEpsilon_E_Float3(orig);
	ray->maxt = maxt;

	ray->time = time;
}

void Ray_Init2_Private(Ray *ray, const float3 orig, const float3 dir, const float time) {
	ray->o.x = orig.x;
	ray->o.y = orig.y;
	ray->o.z = orig.z;

	ray->d.x = dir.x;
	ray->d.y = dir.y;
	ray->d.z = dir.z;

	ray->mint = MachineEpsilon_E_Float3(orig);
	ray->maxt = INFINITY;

	ray->time = time;
}

void Ray_Init4(__global Ray *ray, const float3 orig, const float3 dir,
		const float mint, const float maxt, const float time) {
	VSTORE3F(orig, &ray->o.x);
	VSTORE3F(dir, &ray->d.x);

	ray->mint = mint;
	ray->maxt = maxt;

	ray->time = time;
}

void Ray_Init3(__global Ray *ray, const float3 orig, const float3 dir,
		const float maxt, const float time) {
	VSTORE3F(orig, &ray->o.x);
	VSTORE3F(dir, &ray->d.x);

	ray->mint = MachineEpsilon_E_Float3(orig);
	ray->maxt = maxt;

	ray->time = time;
}

void Ray_Init2(__global Ray *ray, const float3 orig, const float3 dir, const float time) {
	VSTORE3F(orig, &ray->o.x);
	VSTORE3F(dir, &ray->d.x);

	ray->mint = MachineEpsilon_E_Float3(orig);
	ray->maxt = INFINITY;

	ray->time = time;
}

void Ray_ReadAligned4(__global const Ray* restrict ray, float3 *rayOrig, float3 *rayDir,
		float *mint, float *maxt, float *time) {
	__global float4 *basePtr =(__global float4 *)ray;
	const float4 data0 = (*basePtr++);
	const float4 data1 = (*basePtr);

	*rayOrig = (float3)(data0.x, data0.y, data0.z);
	*rayDir = (float3)(data0.w, data1.x, data1.y);

	*mint = data1.z;
	*maxt = data1.w;

	*time = ray->time;
}

void Ray_ReadAligned4_Private(__global const Ray* restrict ray, Ray *dstRay) {
	__global float4 *basePtr =(__global float4 *)ray;
	const float4 data0 = (*basePtr++);
	const float4 data1 = (*basePtr);

	dstRay->o.x = data0.x;
	dstRay->o.y = data0.y;
	dstRay->o.z = data0.z;
	dstRay->d.x = data0.w;
	dstRay->d.y = data1.x;
	dstRay->d.z = data1.y;

	dstRay->mint = data1.z;
	dstRay->maxt = data1.w;

	dstRay->time = ray->time;
}
#line 2 "bbox_types.cl"

/***************************************************************************
 * Copyright 1998-2015 by authors (see AUTHORS.txt)                        *
 *                                                                         *
 *   This file is part of LuxRender.                                       *
 *                                                                         *
 * Licensed under the Apache License, Version 2.0 (the "License");         *
 * you may not use this file except in compliance with the License.        *
 * You may obtain a copy of the License at                                 *
 *                                                                         *
 *     http://www.apache.org/licenses/LICENSE-2.0                          *
 *                                                                         *
 * Unless required by applicable law or agreed to in writing, software     *
 * distributed under the License is distributed on an "AS IS" BASIS,       *
 * WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.*
 * See the License for the specific language governing permissions and     *
 * limitations under the License.                                          *
 ***************************************************************************/

typedef struct {
	Point pMin, pMax;
} BBox;
#line 2 "bbox_funcs.cl"

/***************************************************************************
 * Copyright 1998-2015 by authors (see AUTHORS.txt)                        *
 *                                                                         *
 *   This file is part of LuxRender.                                       *
 *                                                                         *
 * Licensed under the Apache License, Version 2.0 (the "License");         *
 * you may not use this file except in compliance with the License.        *
 * You may obtain a copy of the License at                                 *
 *                                                                         *
 *     http://www.apache.org/licenses/LICENSE-2.0                          *
 *                                                                         *
 * Unless required by applicable law or agreed to in writing, software     *
 * distributed under the License is distributed on an "AS IS" BASIS,       *
 * WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.*
 * See the License for the specific language governing permissions and     *
 * limitations under the License.                                          *
 ***************************************************************************/

int BBox_IntersectP(const float3 pMin, const float3 pMax,
		const float3 rayOrig, const float3 invRayDir,
		const float mint, const float maxt) {
	const float3 l1 = (pMin - rayOrig) * invRayDir;
	const float3 l2 = (pMax - rayOrig) * invRayDir;
	const float3 tNear = fmin(l1, l2);
	const float3 tFar = fmax(l1, l2);

	float t0 = fmax(fmax(fmax(tNear.x, tNear.y), fmax(tNear.x, tNear.z)), mint);
    float t1 = fmin(fmin(fmin(tFar.x, tFar.y), fmin(tFar.x, tFar.z)), maxt);

	return (t1 > t0);
}
#line 2 "transform_types.cl"

/***************************************************************************
 * Copyright 1998-2015 by authors (see AUTHORS.txt)                        *
 *                                                                         *
 *   This file is part of LuxRender.                                       *
 *                                                                         *
 * Licensed under the Apache License, Version 2.0 (the "License");         *
 * you may not use this file except in compliance with the License.        *
 * You may obtain a copy of the License at                                 *
 *                                                                         *
 *     http://www.apache.org/licenses/LICENSE-2.0                          *
 *                                                                         *
 * Unless required by applicable law or agreed to in writing, software     *
 * distributed under the License is distributed on an "AS IS" BASIS,       *
 * WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.*
 * See the License for the specific language governing permissions and     *
 * limitations under the License.                                          *
 ***************************************************************************/

typedef struct {
	Matrix4x4 m, mInv;
} Transform;
#line 2 "transform_funcs.cl"

/***************************************************************************
 * Copyright 1998-2015 by authors (see AUTHORS.txt)                        *
 *                                                                         *
 *   This file is part of LuxRender.                                       *
 *                                                                         *
 * Licensed under the Apache License, Version 2.0 (the "License");         *
 * you may not use this file except in compliance with the License.        *
 * You may obtain a copy of the License at                                 *
 *                                                                         *
 *     http://www.apache.org/licenses/LICENSE-2.0                          *
 *                                                                         *
 * Unless required by applicable law or agreed to in writing, software     *
 * distributed under the License is distributed on an "AS IS" BASIS,       *
 * WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.*
 * See the License for the specific language governing permissions and     *
 * limitations under the License.                                          *
 ***************************************************************************/

float3 Transform_ApplyPoint(__global const Transform* restrict trans, const float3 point) {
	return Matrix4x4_ApplyPoint(&trans->m, point);
}

float3 Transform_ApplyVector(__global const Transform* restrict trans, const float3 vector) {
	return Matrix4x4_ApplyVector(&trans->m, vector);
}

float3 Transform_ApplyNormal(__global const Transform* restrict trans, const float3 normal) {
	return Matrix4x4_ApplyNormal(&trans->m, normal);
}

float3 Transform_InvApplyPoint(__global const Transform* restrict trans, const float3 point) {
	return Matrix4x4_ApplyPoint(&trans->mInv, point);
}

float3 Transform_InvApplyVector(__global const Transform* restrict trans, const float3 vector) {
	return Matrix4x4_ApplyVector(&trans->mInv, vector);
}

float3 Transform_InvApplyNormal(__global const Transform* restrict trans, const float3 normal) {
	return Matrix4x4_ApplyNormal(&trans->mInv, normal);
}
#line 2 "motionsystem_types.cl"

/***************************************************************************
 * Copyright 1998-2015 by authors (see AUTHORS.txt)                        *
 *                                                                         *
 *   This file is part of LuxRender.                                       *
 *                                                                         *
 * Licensed under the Apache License, Version 2.0 (the "License");         *
 * you may not use this file except in compliance with the License.        *
 * You may obtain a copy of the License at                                 *
 *                                                                         *
 *     http://www.apache.org/licenses/LICENSE-2.0                          *
 *                                                                         *
 * Unless required by applicable law or agreed to in writing, software     *
 * distributed under the License is distributed on an "AS IS" BASIS,       *
 * WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.*
 * See the License for the specific language governing permissions and     *
 * limitations under the License.                                          *
 ***************************************************************************/

typedef struct {
	// Scaling
	float Sx, Sy, Sz;
	// Shearing
	float Sxy, Sxz, Syz;
	// Rotation
	Matrix4x4 R;
	// Translation
	float Tx, Ty, Tz;
	// Perspective
	float Px, Py, Pz, Pw;
	// Represents a valid series of transformations
	bool Valid;
} DecomposedTransform;

typedef struct {
	float startTime, endTime;
	Transform start, end;
	DecomposedTransform startT, endT;
	Quaternion startQ, endQ;
	int hasRotation, hasTranslation, hasScale;
	int hasTranslationX, hasTranslationY, hasTranslationZ;
	int hasScaleX, hasScaleY, hasScaleZ;
	// false if start and end transformations are identical
	int isActive;
} InterpolatedTransform;

typedef struct {
	unsigned int interpolatedTransformFirstIndex;
	unsigned int interpolatedTransformLastIndex;
} MotionSystem;
#line 2 "motionsystem_funcs.cl"

/***************************************************************************
 * Copyright 1998-2015 by authors (see AUTHORS.txt)                        *
 *                                                                         *
 *   This file is part of LuxRender.                                       *
 *                                                                         *
 * Licensed under the Apache License, Version 2.0 (the "License");         *
 * you may not use this file except in compliance with the License.        *
 * You may obtain a copy of the License at                                 *
 *                                                                         *
 *     http://www.apache.org/licenses/LICENSE-2.0                          *
 *                                                                         *
 * Unless required by applicable law or agreed to in writing, software     *
 * distributed under the License is distributed on an "AS IS" BASIS,       *
 * WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.*
 * See the License for the specific language governing permissions and     *
 * limitations under the License.                                          *
 ***************************************************************************/

void InterpolatedTransform_Sample(__global const InterpolatedTransform* restrict interpolatedTransform,
		const float time, Matrix4x4 *result) {
	if (!interpolatedTransform->isActive) {
		*result = interpolatedTransform->start.m;
		return;
	}

	// Determine interpolation value
	if (time <= interpolatedTransform->startTime) {
		*result = interpolatedTransform->start.m;
		return;
	}
	if (time >= interpolatedTransform->endTime) {
		*result = interpolatedTransform->end.m;
		return;
	}

	const float w = interpolatedTransform->endTime - interpolatedTransform->startTime;
	const float d = time - interpolatedTransform->startTime;
	const float le = d / w;

	// if translation only, just modify start matrix
	if (interpolatedTransform->hasTranslation &&
			!(interpolatedTransform->hasScale || interpolatedTransform->hasRotation)) {
		*result = interpolatedTransform->start.m;
		if (interpolatedTransform->hasTranslationX)
			result->m[0][3] = mix(interpolatedTransform->startT.Tx, interpolatedTransform->endT.Tx, le);
		if (interpolatedTransform->hasTranslationY)
			result->m[1][3] = mix(interpolatedTransform->startT.Ty, interpolatedTransform->endT.Ty, le);
		if (interpolatedTransform->hasTranslationZ)
			result->m[2][3] = mix(interpolatedTransform->startT.Tz, interpolatedTransform->endT.Tz, le);

		return;
	}

	if (interpolatedTransform->hasRotation) {
		// Quaternion interpolation of rotation
		const float4 startQ = VLOAD4F(&interpolatedTransform->startQ.w);
		const float4 endQ = VLOAD4F(&interpolatedTransform->endQ.w);
		const float4 interQ = Quaternion_Slerp(le, startQ, endQ);
		Quaternion_ToMatrix(interQ, result);
	} else
		*result = interpolatedTransform->startT.R;

	if (interpolatedTransform->hasScale) {
		const float Sx = mix(interpolatedTransform->startT.Sx, interpolatedTransform->endT.Sx, le);
		const float Sy = mix(interpolatedTransform->startT.Sy, interpolatedTransform->endT.Sy, le); 
		const float Sz = mix(interpolatedTransform->startT.Sz, interpolatedTransform->endT.Sz, le);

		// T * S * R
		for (uint j = 0; j < 3; ++j) {
			result->m[0][j] = Sx * result->m[0][j];
			result->m[1][j] = Sy * result->m[1][j];
			result->m[2][j] = Sz * result->m[2][j];
		}
	} else {
		for (uint j = 0; j < 3; ++j) {
			result->m[0][j] = interpolatedTransform->startT.Sx * result->m[0][j];
			result->m[1][j] = interpolatedTransform->startT.Sy * result->m[1][j];
			result->m[2][j] = interpolatedTransform->startT.Sz * result->m[2][j];
		}
	}

	if (interpolatedTransform->hasTranslationX)
		result->m[0][3] = mix(interpolatedTransform->startT.Tx, interpolatedTransform->endT.Tx, le);
	else
		result->m[0][3] = interpolatedTransform->startT.Tx;

	if (interpolatedTransform->hasTranslationY)
		result->m[1][3] = mix(interpolatedTransform->startT.Ty, interpolatedTransform->endT.Ty, le);
	else
		result->m[1][3] = interpolatedTransform->startT.Ty;

	if (interpolatedTransform->hasTranslationZ)
		result->m[2][3] = mix(interpolatedTransform->startT.Tz, interpolatedTransform->endT.Tz, le);
	else
		result->m[2][3] = interpolatedTransform->startT.Tz;
}

void MotionSystem_Sample(__global const MotionSystem* restrict motionSystem, const float time,
		__global const InterpolatedTransform *interpolatedTransforms, Matrix4x4 *result) {
	const uint interpolatedTransformFirstIndex = motionSystem->interpolatedTransformFirstIndex;
	const uint interpolatedTransformLastIndex = motionSystem->interpolatedTransformLastIndex;

	// Pick the right InterpolatedTransform
	uint index = interpolatedTransformLastIndex;
	for (uint i = interpolatedTransformFirstIndex; i <= interpolatedTransformLastIndex; ++i) {
		if (time < interpolatedTransforms[i].endTime) {
			index = i;
			break;
		}
	}

	InterpolatedTransform_Sample(&interpolatedTransforms[index], time, result);
}
#line 2 "qbvh_types.cl"

/***************************************************************************
 * Copyright 1998-2015 by authors (see AUTHORS.txt)                        *
 *                                                                         *
 *   This file is part of LuxRender.                                       *
 *                                                                         *
 * Licensed under the Apache License, Version 2.0 (the "License");         *
 * you may not use this file except in compliance with the License.        *
 * You may obtain a copy of the License at                                 *
 *                                                                         *
 *     http://www.apache.org/licenses/LICENSE-2.0                          *
 *                                                                         *
 * Unless required by applicable law or agreed to in writing, software     *
 * distributed under the License is distributed on an "AS IS" BASIS,       *
 * WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.*
 * See the License for the specific language governing permissions and     *
 * limitations under the License.                                          *
 ***************************************************************************/

typedef struct QuadRay {
	float4 ox, oy, oz;
	float4 dx, dy, dz;
	float4 mint, maxt;
} QuadRay;

typedef struct {
	float4 origx, origy, origz;
	float4 edge1x, edge1y, edge1z;
	float4 edge2x, edge2y, edge2z;
	uint4 meshIndex, triangleIndex;
} QuadTiangle;

typedef struct {
	float4 bboxes[2][3];
	int4 children;
} QBVHNode;

#define emptyLeafNode 0xffffffff

#define QBVHNode_IsLeaf(index) (index < 0)
#define QBVHNode_IsEmpty(index) (index == emptyLeafNode)
#define QBVHNode_NbQuadPrimitives(index) ((uint)(((index >> 27) & 0xf) + 1))
#define QBVHNode_FirstQuadIndex(index) (index & 0x07ffffff)

// Using invDir0/invDir1/invDir2 instead of an
// array because I dont' trust OpenCL compiler =)
int4 QBVHNode_BBoxIntersect(
        const float4 bboxes_minX, const float4 bboxes_maxX,
        const float4 bboxes_minY, const float4 bboxes_maxY,
        const float4 bboxes_minZ, const float4 bboxes_maxZ,
        const QuadRay *ray4,
		const float4 invDir0, const float4 invDir1, const float4 invDir2) {
	float4 tMin = ray4->mint;
	float4 tMax = ray4->maxt;

	// X coordinate
	tMin = fmax(tMin, (bboxes_minX - ray4->ox) * invDir0);
	tMax = fmin(tMax, (bboxes_maxX - ray4->ox) * invDir0);

	// Y coordinate
	tMin = fmax(tMin, (bboxes_minY - ray4->oy) * invDir1);
	tMax = fmin(tMax, (bboxes_maxY - ray4->oy) * invDir1);

	// Z coordinate
	tMin = fmax(tMin, (bboxes_minZ - ray4->oz) * invDir2);
	tMax = fmin(tMax, (bboxes_maxZ - ray4->oz) * invDir2);

	// Return the visit flags
	return  (tMax >= tMin);
}

void QuadTriangle_Intersect(
    const float4 origx, const float4 origy, const float4 origz,
    const float4 edge1x, const float4 edge1y, const float4 edge1z,
    const float4 edge2x, const float4 edge2y, const float4 edge2z,
    const uint4 meshIndex,  const uint4 triangleIndex,
    QuadRay *ray4, RayHit *rayHit) {
	//--------------------------------------------------------------------------
	// Calc. b1 coordinate

	const float4 s1x = (ray4->dy * edge2z) - (ray4->dz * edge2y);
	const float4 s1y = (ray4->dz * edge2x) - (ray4->dx * edge2z);
	const float4 s1z = (ray4->dx * edge2y) - (ray4->dy * edge2x);

	const float4 divisor = (s1x * edge1x) + (s1y * edge1y) + (s1z * edge1z);

	const float4 dx = ray4->ox - origx;
	const float4 dy = ray4->oy - origy;
	const float4 dz = ray4->oz - origz;

	const float4 b1 = ((dx * s1x) + (dy * s1y) + (dz * s1z)) / divisor;

	//--------------------------------------------------------------------------
	// Calc. b2 coordinate

	const float4 s2x = (dy * edge1z) - (dz * edge1y);
	const float4 s2y = (dz * edge1x) - (dx * edge1z);
	const float4 s2z = (dx * edge1y) - (dy * edge1x);

	const float4 b2 = ((ray4->dx * s2x) + (ray4->dy * s2y) + (ray4->dz * s2z)) / divisor;

	//--------------------------------------------------------------------------
	// Calc. b0 coordinate

	const float4 b0 = ((float4)1.f) - b1 - b2;

	//--------------------------------------------------------------------------

	const float4 t = ((edge2x * s2x) + (edge2y * s2y) + (edge2z * s2z)) / divisor;

    float _b1, _b2;
	float maxt = ray4->maxt.s0;
    uint mIndex, tIndex;

    int4 cond = isnotequal(divisor, (float4)0.f) & isgreaterequal(b0, (float4)0.f) &
			isgreaterequal(b1, (float4)0.f) & isgreaterequal(b2, (float4)0.f) &
			isgreater(t, ray4->mint);

    const int cond0 = cond.s0 && (t.s0 < maxt);
    maxt = select(maxt, t.s0, cond0);
    _b1 = select(0.f, b1.s0, cond0);
    _b2 = select(0.f, b2.s0, cond0);
    mIndex = select(NULL_INDEX, meshIndex.s0, cond0);
	tIndex = select(NULL_INDEX, triangleIndex.s0, cond0);

    const int cond1 = cond.s1 && (t.s1 < maxt);
    maxt = select(maxt, t.s1, cond1);
    _b1 = select(_b1, b1.s1, cond1);
    _b2 = select(_b2, b2.s1, cond1);
    mIndex = select(mIndex, meshIndex.s1, cond1);
	tIndex = select(tIndex, triangleIndex.s1, cond1);

    const int cond2 = cond.s2 && (t.s2 < maxt);
    maxt = select(maxt, t.s2, cond2);
    _b1 = select(_b1, b1.s2, cond2);
    _b2 = select(_b2, b2.s2, cond2);
    mIndex = select(mIndex, meshIndex.s2, cond2);
	tIndex = select(tIndex, triangleIndex.s2, cond2);

    const int cond3 = cond.s3 && (t.s3 < maxt);
    maxt = select(maxt, t.s3, cond3);
    _b1 = select(_b1, b1.s3, cond3);
    _b2 = select(_b2, b2.s3, cond3);
    mIndex = select(mIndex, meshIndex.s3, cond3);
	tIndex = select(tIndex, triangleIndex.s3, cond3);

	if (mIndex == NULL_INDEX)
		return;

	ray4->maxt = (float4)maxt;

	rayHit->t = maxt;
	rayHit->b1 = _b1;
	rayHit->b2 = _b2;
	rayHit->meshIndex = mIndex;
	rayHit->triangleIndex = tIndex;
}
#line 2 "mqbvh_kernel.cl"

/***************************************************************************
 * Copyright 1998-2015 by authors (see AUTHORS.txt)                        *
 *                                                                         *
 *   This file is part of LuxRender.                                       *
 *                                                                         *
 * Licensed under the Apache License, Version 2.0 (the "License");         *
 * you may not use this file except in compliance with the License.        *
 * You may obtain a copy of the License at                                 *
 *                                                                         *
 *     http://www.apache.org/licenses/LICENSE-2.0                          *
 *                                                                         *
 * Unless required by applicable law or agreed to in writing, software     *
 * distributed under the License is distributed on an "AS IS" BASIS,       *
 * WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.*
 * See the License for the specific language governing permissions and     *
 * limitations under the License.                                          *
 ***************************************************************************/

// Using a large stack size to avoid the allocation of the array on
// GPU registers (otherwise the GPU can easily run out of registers)
#define STACK_SIZE 64

void LeafIntersect(
		const Ray *ray,
		RayHit *rayHit,
		__global const QBVHNode* restrict nodes,
		__global const QuadTiangle* restrict quadTris) {
	// Prepare the ray for intersection
	QuadRay ray4;
    ray4.ox = (float4)ray->o.x;
    ray4.oy = (float4)ray->o.y;
    ray4.oz = (float4)ray->o.z;

    ray4.dx = (float4)ray->d.x;
    ray4.dy = (float4)ray->d.y;
    ray4.dz = (float4)ray->d.z;

    ray4.mint = (float4)ray->mint;
    ray4.maxt = (float4)ray->maxt;

	const float4 invDir0 = (float4)(1.f / ray4.dx.s0);
	const float4 invDir1 = (float4)(1.f / ray4.dy.s0);
	const float4 invDir2 = (float4)(1.f / ray4.dz.s0);

	const int signs0 = signbit(ray4.dx.s0);
	const int signs1 = signbit(ray4.dy.s0);
	const int signs2 = signbit(ray4.dz.s0);

	rayHit->meshIndex = NULL_INDEX;
	rayHit->triangleIndex = NULL_INDEX;

	//------------------------------
	// Main loop
	int todoNode = 0; // the index in the stack
	int nodeStack[STACK_SIZE];
	nodeStack[0] = 0; // first node to handle: root node

	while (todoNode >= 0) {
		const int nodeData = nodeStack[todoNode];
		--todoNode;

		// Leaves are identified by a negative index
		if (!QBVHNode_IsLeaf(nodeData)) {
			__global const QBVHNode* restrict node = &nodes[nodeData];
            const int4 visit = QBVHNode_BBoxIntersect(
                node->bboxes[signs0][0], node->bboxes[1 - signs0][0],
                node->bboxes[signs1][1], node->bboxes[1 - signs1][1],
                node->bboxes[signs2][2], node->bboxes[1 - signs2][2],
                &ray4,
				invDir0, invDir1, invDir2);

			const int4 children = node->children;

			// For some reason doing logic operations with int4 is very slow
			nodeStack[todoNode + 1] = children.s3;
			todoNode += (visit.s3 && !QBVHNode_IsEmpty(children.s3)) ? 1 : 0;
			nodeStack[todoNode + 1] = children.s2;
			todoNode += (visit.s2 && !QBVHNode_IsEmpty(children.s2)) ? 1 : 0;
			nodeStack[todoNode + 1] = children.s1;
			todoNode += (visit.s1 && !QBVHNode_IsEmpty(children.s1)) ? 1 : 0;
			nodeStack[todoNode + 1] = children.s0;
			todoNode += (visit.s0 && !QBVHNode_IsEmpty(children.s0)) ? 1 : 0;
		} else {
			// Perform intersection
			const uint nbQuadPrimitives = QBVHNode_NbQuadPrimitives(nodeData);
			const uint offset = QBVHNode_FirstQuadIndex(nodeData);

			for (uint primNumber = offset; primNumber < (offset + nbQuadPrimitives); ++primNumber) {
                __global const QuadTiangle* restrict quadTri = &quadTris[primNumber];
                const float4 origx = quadTri->origx;
                const float4 origy = quadTri->origy;
                const float4 origz = quadTri->origz;
                const float4 edge1x = quadTri->edge1x;
                const float4 edge1y = quadTri->edge1y;
                const float4 edge1z = quadTri->edge1z;
                const float4 edge2x = quadTri->edge2x;
                const float4 edge2y = quadTri->edge2y;
                const float4 edge2z = quadTri->edge2z;
                const uint4 meshIndex = quadTri->meshIndex;
				const uint4 triangleIndex = quadTri->triangleIndex;

				QuadTriangle_Intersect(
                    origx, origy, origz,
                    edge1x, edge1y, edge1z,
                    edge2x, edge2y, edge2z,
                    meshIndex, triangleIndex,
                    &ray4, rayHit);
            }
		}
	}
}

#define MQBVH_TRANSFORMATIONS_PARAM_DECL , __global const uint* restrict leafTransformationIndex , __global const Matrix4x4* restrict leafTransformations
#define MQBVH_TRANSFORMATIONS_PARAM , leafTransformationIndex, leafTransformations

#define MQBVH_MOTIONSYSTEMS_PARAM_DECL , __global const MotionSystem* restrict leafMotionSystems , __global const InterpolatedTransform* restrict leafInterpolatedTransforms
#define MQBVH_MOTIONSYSTEMS_PARAM , leafMotionSystems, leafInterpolatedTransforms

#define ACCELERATOR_INTERSECT_PARAM_DECL ,__global const QBVHNode* restrict nodes, __global const uint* restrict qbvhMemMap, __global const QBVHNode* restrict leafNodes, __global const QuadTiangle* restrict const leafQuadTris MQBVH_TRANSFORMATIONS_PARAM_DECL MQBVH_MOTIONSYSTEMS_PARAM_DECL
#define ACCELERATOR_INTERSECT_PARAM ,nodes, qbvhMemMap, leafNodes, leafQuadTris MQBVH_TRANSFORMATIONS_PARAM MQBVH_MOTIONSYSTEMS_PARAM

void Accelerator_Intersect(
		const Ray *ray,
		RayHit *rayHit
		ACCELERATOR_INTERSECT_PARAM_DECL
		) {
	// Prepare the ray for intersection
    const float3 rayOrig = VLOAD3F_Private(&ray->o.x);
    const float3 rayDir = VLOAD3F_Private(&ray->d.x);

	QuadRay ray4;
	ray4.ox = (float4)ray->o.x;
	ray4.oy = (float4)ray->o.y;
	ray4.oz = (float4)ray->o.z;

	ray4.dx = (float4)ray->d.x;
	ray4.dy = (float4)ray->d.y;
	ray4.dz = (float4)ray->d.z;

	ray4.mint = (float4)ray->mint;
	ray4.maxt = (float4)ray->maxt;

	const float4 invDir0 = (float4)(1.f / ray4.dx.s0);
	const float4 invDir1 = (float4)(1.f / ray4.dy.s0);
	const float4 invDir2 = (float4)(1.f / ray4.dz.s0);

	const int signs0 = signbit(ray4.dx.s0);
	const int signs1 = signbit(ray4.dy.s0);
	const int signs2 = signbit(ray4.dz.s0);

	const float rayTime = ray->time;

	rayHit->meshIndex = NULL_INDEX;
	rayHit->triangleIndex = NULL_INDEX;

	//------------------------------
	// Main loop
	int todoNode = 0; // the index in the stack
	int nodeStack[STACK_SIZE];
	nodeStack[0] = 0; // first node to handle: root node

	while (todoNode >= 0) {
		const int nodeData = nodeStack[todoNode];
		--todoNode;

		// Leaves are identified by a negative index
		if (!QBVHNode_IsLeaf(nodeData)) {
			__global const QBVHNode* restrict node = &nodes[nodeData];
            const int4 visit = QBVHNode_BBoxIntersect(
                node->bboxes[signs0][0], node->bboxes[1 - signs0][0],
                node->bboxes[signs1][1], node->bboxes[1 - signs1][1],
                node->bboxes[signs2][2], node->bboxes[1 - signs2][2],
                &ray4,
				invDir0, invDir1, invDir2);

			const int4 children = node->children;

			// For some reason doing logic operations with int4 are very slow
			nodeStack[todoNode + 1] = children.s3;
			todoNode += (visit.s3 && !QBVHNode_IsEmpty(children.s3)) ? 1 : 0;
			nodeStack[todoNode + 1] = children.s2;
			todoNode += (visit.s2 && !QBVHNode_IsEmpty(children.s2)) ? 1 : 0;
			nodeStack[todoNode + 1] = children.s1;
			todoNode += (visit.s1 && !QBVHNode_IsEmpty(children.s1)) ? 1 : 0;
			nodeStack[todoNode + 1] = children.s0;
			todoNode += (visit.s0 && !QBVHNode_IsEmpty(children.s0)) ? 1 : 0;
		} else {
			// Perform intersection with QBVH leaf
			const uint leafIndex = QBVHNode_FirstQuadIndex(nodeData);

            Ray tray;
			float3 newRayOrig = rayOrig;
			float3 newRayDir = rayDir;

			if (leafTransformationIndex[leafIndex] != NULL_INDEX) {
				__global const Matrix4x4* restrict m = &leafTransformations[leafTransformationIndex[leafIndex]];
				newRayOrig = Matrix4x4_ApplyPoint_Align(m, newRayOrig);
				newRayDir = Matrix4x4_ApplyVector(m, newRayDir);
			}

			if (leafMotionSystems[leafIndex].interpolatedTransformFirstIndex != NULL_INDEX) {
				// Transform ray origin and direction
				Matrix4x4 m;
				MotionSystem_Sample(&leafMotionSystems[leafIndex], rayTime, leafInterpolatedTransforms, &m);
				newRayOrig = Matrix4x4_ApplyPoint_Private(&m, newRayOrig);
				newRayDir = Matrix4x4_ApplyVector_Private(&m, newRayDir);
			}

			tray.o.x = newRayOrig.x;
			tray.o.y = newRayOrig.y;
			tray.o.z = newRayOrig.z;
			tray.d.x = newRayDir.x;
			tray.d.y = newRayDir.y;
			tray.d.z = newRayDir.z;
			tray.mint = ray4.mint.s0;
			tray.maxt = ray4.maxt.s0;
			tray.time = rayTime;

            const uint memIndex = leafIndex * 2;
            const uint leafNodeOffset = qbvhMemMap[memIndex];
            __global const QBVHNode* restrict n = &leafNodes[leafNodeOffset];
            const uint leafQuadTriOffset = qbvhMemMap[memIndex + 1];
            __global const QuadTiangle* restrict qt = &leafQuadTris[leafQuadTriOffset];

            RayHit tmpRayHit;
            LeafIntersect(&tray, &tmpRayHit, n, qt);

            if (tmpRayHit.meshIndex != NULL_INDEX) {
                rayHit->t = tmpRayHit.t;
                rayHit->b1 = tmpRayHit.b1;
                rayHit->b2 = tmpRayHit.b2;
                rayHit->meshIndex = leafIndex;
				rayHit->triangleIndex = tmpRayHit.triangleIndex;

                ray4.maxt = (float4)tmpRayHit.t;
            }
		}
	}
}

__kernel __attribute__((work_group_size_hint(64, 1, 1))) void Accelerator_Intersect_RayBuffer(
		__global const Ray* restrict rays,
		__global RayHit *rayHits,
		const uint rayCount
		ACCELERATOR_INTERSECT_PARAM_DECL
		) {
	// Select the ray to check
	const int gid = get_global_id(0);
	if (gid >= rayCount)
		return;

	Ray ray;
	Ray_ReadAligned4_Private(&rays[gid], &ray);

	RayHit rayHit;
	Accelerator_Intersect(
		&ray,
		&rayHit
		ACCELERATOR_INTERSECT_PARAM
		);

	// Write result
	__global RayHit *memRayHit = &rayHits[gid];
	memRayHit->t = rayHit.t;
	memRayHit->b1 = rayHit.b1;
	memRayHit->b2 = rayHit.b2;
	memRayHit->meshIndex = rayHit.meshIndex;
	memRayHit->triangleIndex = rayHit.triangleIndex;
}
