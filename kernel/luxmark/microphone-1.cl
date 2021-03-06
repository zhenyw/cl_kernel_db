#define QBVH_STACK_SIZE 48
#define USE_IMAGE_STORAGE
#define QBVH_USE_LOCAL_MEMORY
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
#line 2 "qbvh_kernel.cl"

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

#if defined(QBVH_USE_LOCAL_MEMORY)
#define QBVH_LOCAL_MEMORY_PARAM_DECL , __local int *nodeStacks
#define QBVH_LOCAL_MEMORY_PARAM , nodeStacks
#else
#define QBVH_LOCAL_MEMORY_PARAM_DECL
#define QBVH_LOCAL_MEMORY_PARAM
#endif

#ifdef USE_IMAGE_STORAGE
#define ACCELERATOR_INTERSECT_PARAM_DECL , __read_only image2d_t nodes, __read_only image2d_t quadTris QBVH_LOCAL_MEMORY_PARAM_DECL
#define ACCELERATOR_INTERSECT_PARAM , nodes, quadTris QBVH_LOCAL_MEMORY_PARAM
#else
#define ACCELERATOR_INTERSECT_PARAM_DECL ,__global const QBVHNode* restrict nodes, __global const QuadTiangle* restrict quadTris QBVH_LOCAL_MEMORY_PARAM_DECL
#define ACCELERATOR_INTERSECT_PARAM , nodes, quadTris QBVH_LOCAL_MEMORY_PARAM
#endif

void Accelerator_Intersect(
		const Ray *ray,
		RayHit *rayHit
		ACCELERATOR_INTERSECT_PARAM_DECL
		) {
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

	const int isigns0 = 1 - signs0;
	const int isigns1 = 1 - signs1;
	const int isigns2 = 1 - signs2;

	rayHit->meshIndex = NULL_INDEX;
	rayHit->triangleIndex = NULL_INDEX;

	//------------------------------
	// Main loop
	int todoNode = 0; // the index in the stack
	// nodeStack leads to a lot of local memory banks conflicts however it has not real
	// impact on performances (I guess access latency is hidden by other stuff).
	// Avoiding conflicts is easy to do but it requires to know the work group
	// size (not worth doing if there are not performance benefits).
#if defined(QBVH_USE_LOCAL_MEMORY)
	__local int *nodeStack = &nodeStacks[QBVH_STACK_SIZE * get_local_id(0)];
#else
	int nodeStack[QBVH_STACK_SIZE];
#endif
	nodeStack[0] = 0; // first node to handle: root node

#ifdef USE_IMAGE_STORAGE
    const int quadTrisImageWidth = get_image_width(quadTris);

    const int bboxes_minXIndex = (signs0 * 3);
    const int bboxes_maxXIndex = (isigns0 * 3);
    const int bboxes_minYIndex = (signs1 * 3) + 1;
    const int bboxes_maxYIndex = (isigns1 * 3) + 1;
    const int bboxes_minZIndex = (signs2 * 3) + 2;
    const int bboxes_maxZIndex = (isigns2 * 3) + 2;

    const sampler_t imageSampler = CLK_NORMALIZED_COORDS_FALSE | CLK_ADDRESS_CLAMP | CLK_FILTER_NEAREST;
#endif

	//int maxDepth = 0;
	while (todoNode >= 0) {
		const int nodeData = nodeStack[todoNode];
		--todoNode;

		// Leaves are identified by a negative index
		if (!QBVHNode_IsLeaf(nodeData)) {
#ifdef USE_IMAGE_STORAGE
            // Read the node information from the image storage

			// 7 pixels required for the storage of a QBVH node
            const ushort inx = (nodeData >> 16) * 7;
            const ushort iny = (nodeData & 0xffff);
            const float4 bboxes_minX = as_float4(read_imageui(nodes, imageSampler, (int2)(inx + bboxes_minXIndex, iny)));
            const float4 bboxes_maxX = as_float4(read_imageui(nodes, imageSampler, (int2)(inx + bboxes_maxXIndex, iny)));
            const float4 bboxes_minY = as_float4(read_imageui(nodes, imageSampler, (int2)(inx + bboxes_minYIndex, iny)));
            const float4 bboxes_maxY = as_float4(read_imageui(nodes, imageSampler, (int2)(inx + bboxes_maxYIndex, iny)));
            const float4 bboxes_minZ = as_float4(read_imageui(nodes, imageSampler, (int2)(inx + bboxes_minZIndex, iny)));
            const float4 bboxes_maxZ = as_float4(read_imageui(nodes, imageSampler, (int2)(inx + bboxes_maxZIndex, iny)));
            const int4 children = as_int4(read_imageui(nodes, imageSampler, (int2)(inx + 6, iny)));

			const int4 visit = QBVHNode_BBoxIntersect(
                bboxes_minX, bboxes_maxX,
                bboxes_minY, bboxes_maxY,
                bboxes_minZ, bboxes_maxZ,
                &ray4,
				invDir0, invDir1, invDir2);
#else
			__global const QBVHNode* restrict node = &nodes[nodeData];
            const int4 visit = QBVHNode_BBoxIntersect(
                node->bboxes[signs0][0], node->bboxes[isigns0][0],
                node->bboxes[signs1][1], node->bboxes[isigns1][1],
                node->bboxes[signs2][2], node->bboxes[isigns2][2],
                &ray4,
				invDir0, invDir1, invDir2);

			const int4 children = node->children;
#endif

			// For some reason doing logic operations with int4 is very slow
			nodeStack[todoNode + 1] = children.s3;
			todoNode += (visit.s3 && !QBVHNode_IsEmpty(children.s3)) ? 1 : 0;
			nodeStack[todoNode + 1] = children.s2;
			todoNode += (visit.s2 && !QBVHNode_IsEmpty(children.s2)) ? 1 : 0;
			nodeStack[todoNode + 1] = children.s1;
			todoNode += (visit.s1 && !QBVHNode_IsEmpty(children.s1)) ? 1 : 0;
			nodeStack[todoNode + 1] = children.s0;
			todoNode += (visit.s0 && !QBVHNode_IsEmpty(children.s0)) ? 1 : 0;

			//maxDepth = max(maxDepth, todoNode);
		} else {
			// Perform intersection
			const uint nbQuadPrimitives = QBVHNode_NbQuadPrimitives(nodeData);
			const uint offset = QBVHNode_FirstQuadIndex(nodeData);

#ifdef USE_IMAGE_STORAGE
			// 11 pixels required for the storage of QBVH Triangles
            ushort inx = (offset >> 16) * 11;
            ushort iny = (offset & 0xffff);
#endif

			for (uint primNumber = offset; primNumber < (offset + nbQuadPrimitives); ++primNumber) {
#ifdef USE_IMAGE_STORAGE
                const float4 origx = as_float4(read_imageui(quadTris, imageSampler, (int2)(inx++, iny)));
                const float4 origy = as_float4(read_imageui(quadTris, imageSampler, (int2)(inx++, iny)));
                const float4 origz = as_float4(read_imageui(quadTris, imageSampler, (int2)(inx++, iny)));
                const float4 edge1x = as_float4(read_imageui(quadTris, imageSampler, (int2)(inx++, iny)));
                const float4 edge1y = as_float4(read_imageui(quadTris, imageSampler, (int2)(inx++, iny)));
                const float4 edge1z = as_float4(read_imageui(quadTris, imageSampler, (int2)(inx++, iny)));
                const float4 edge2x = as_float4(read_imageui(quadTris, imageSampler, (int2)(inx++, iny)));
                const float4 edge2y = as_float4(read_imageui(quadTris, imageSampler, (int2)(inx++, iny)));
                const float4 edge2z = as_float4(read_imageui(quadTris, imageSampler, (int2)(inx++, iny)));
                const uint4 meshIndex = read_imageui(quadTris, imageSampler, (int2)(inx++, iny));
				const uint4 triangleIndex = read_imageui(quadTris, imageSampler, (int2)(inx++, iny));

                if (inx >= quadTrisImageWidth) {
                    inx = 0;
                    iny++;
                }
#else
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
#endif
				QuadTriangle_Intersect(
                    origx, origy, origz,
                    edge1x, edge1y, edge1z,
                    edge2x, edge2y, edge2z,
                    meshIndex, triangleIndex,
                    &ray4, rayHit);
            }
		}
	}

	//printf("MaxDepth=%02d\n", maxDepth);
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
