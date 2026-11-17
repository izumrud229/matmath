# matmath

**Header only library for Vulkan and OpenGL API on C11**

## Features
- **Header-only**: Just drop `matmath.h` into your project and include it.
- **GPU-ready memory alignment**: `vec4`, `m4`, and `quat` are 16-byte aligned (`_Alignas(16)`) for direct `memcpy` transfers into UBOs, SSBOs, and Push Constants.
- **Vulkan & OpenGL support**: Dedicated projection functions (`m4_perspective_vk`, `m4_ortho_vk`) supporting Y-axis flipping (`flip_y`) and `[0, 1]` depth range.
- **Zero external dependencies**: Relies strictly on the C standard library (`<math.h>`, `<stdalign.h>`).
- **Weight**: 18KB header-only.

## What's Included
- **Scalars**: `clampf`, `lerpf`, `radiansf`, `degreesf`.
- **Vectors**: `vec2`, `vec3`, `vec4` (addition, subtraction, dot/cross products, normalization, linear interpolation).
- **Matrices**: `m2`, `m3`, `m4` (transformations, multiplication, transposition, inversion, normal matrix `m3_normal_matrix`, view matrix `m4_look_at`).
- **Quaternions**: `quat` (multiplication, normalization, spherical linear interpolation `qslerp`, matrix conversion `q_to_m4`).

## Quick Start

```c
#include "matmath.h"

int main(void) {
    vec3 eye = vmake(0.0f, 2.0f, 5.0f);
    vec3 center = vmake(0.0f, 0.0f, 0.0f);
    vec3 up = vmake(0.0f, 1.0f, 0.0f);

    m4 view;
    m4_look_at(&view, eye, center, up);

    m4 proj;
    m4_perspective_vk(&proj, radiansf(60.0f), 16.0f / 9.0f, 0.1f, 100.0f, 1);

    return 0;
}
