#ifndef MATMATH_H
#define MATMATH_H

#include <math.h>
#include <stdalign.h>

//--++=== structures ===++--

typedef struct { float x, y; } vec2;
typedef struct { float x, y, z; } vec3;
typedef struct { _Alignas(16) float x, y, z, w; } vec4;
typedef struct { float m[4]; } m2;
typedef struct { float m[9]; } m3;
typedef struct { _Alignas(16) float m[16]; } m4;
typedef struct { _Alignas(16) float x, y, z, w; } quat;

//--++=== scalar ==++--

#define MLIB_PI 3.14159265358979323846f

static inline float clampf(float v, float lo, float hi) {
    return v < lo ? lo : (v > hi ? hi : v);
}

static inline float lerpf(float a, float b, float t) {
    return a + t * (b - a);
}

static inline float radiansf(float deg) {
    return deg * (MLIB_PI / 180.0f);
}

static inline float degreesf(float rad) {
    return rad * (180.0f / MLIB_PI);
}

//--++=== vectors(x, y, z) ==++--

static inline vec3 vmake(float x, float y, float z) {
    vec3 v = {x, y, z};
    return v;
}

static inline vec3 vadd(vec3 a, vec3 b) {
    return vmake(a.x + b.x, a.y + b.y, a.z + b.z);
}

static inline vec3 vsub(vec3 a, vec3 b) {
    return vmake(a.x - b.x, a.y - b.y, a.z - b.z);
}

static inline vec3 vscale(vec3 v, float s) {
    return vmake(v.x * s, v.y * s, v.z * s);
}

static inline float vdot(vec3 a, vec3 b) {
    return a.x * b.x + a.y * b.y + a.z * b.z;
}

static inline vec3 vcross(vec3 a, vec3 b) {
    return vmake(
        a.y * b.z - a.z * b.y,
        a.z * b.x - a.x * b.z,
        a.x * b.y - a.y * b.x
    );
}

static inline float vlength_sq(vec3 v) {
    return vdot(v, v);
}

static inline float vlength(vec3 v) {
    return sqrtf(vlength_sq(v));
}

static inline vec3 vnormalize(vec3 v) {
    float len = vlength(v);
    if (len < 1e-8f) return vmake(0, 0, 0);
    return vscale(v, 1.0f / len);
}

static inline vec3 vlerp(vec3 a, vec3 b, float t) {
    return vmake(a.x + t*(b.x-a.x), a.y + t*(b.y-a.y), a.z + t*(b.z-a.z));
}

//--++=== matrix 4x4 ===++--

static inline void m4_identity(m4 *out) {
    for (int i = 0; i < 16; i++) out->m[i] = 0.0f;
    out->m[0] = out->m[5] = out->m[10] = out->m[15] = 1.0f;
}

static inline void m4_translate(m4 *out, float tx, float ty, float tz) {
    m4_identity(out);
    out->m[12] = tx;
    out->m[13] = ty;
    out->m[14] = tz;
}

static inline void m4_scale(m4 *out, float sx, float sy, float sz) {
    m4_identity(out);
    out->m[0]  = sx;
    out->m[5]  = sy;
    out->m[10] = sz;
}

static inline void m4_rotate_x(m4 *out, float angle) {
    float c = cosf(angle), s = sinf(angle);
    m4_identity(out);
    out->m[5]  =  c;
    out->m[6]  =  s;
    out->m[9]  = -s;
    out->m[10] =  c;
}

static inline void m4_rotate_y(m4 *out, float angle) {
    float c = cosf(angle), s = sinf(angle);
    m4_identity(out);
    out->m[0]  =  c;
    out->m[2]  = -s;
    out->m[8]  =  s;
    out->m[10] =  c;
}

static inline void m4_rotate_z(m4 *out, float angle) {
    float c = cosf(angle), s = sinf(angle);
    m4_identity(out);
    out->m[0] =  c;
    out->m[1] =  s;
    out->m[4] = -s;
    out->m[5] =  c;
}

static inline void m4_multiply(m4 *out, const m4 *a, const m4 *b) {
    m4 tmp;
    for (int col = 0; col < 4; col++) {
        for (int row = 0; row < 4; row++) {
            float sum = 0.0f;
            for (int k = 0; k < 4; k++)
                sum += a->m[k * 4 + row] * b->m[col * 4 + k];
            tmp.m[col * 4 + row] = sum;
        }
    }
    *out = tmp;
}

static inline void m4_mul_vec4(float out[4], const m4 *m, const float v[4]) {
    float tmp[4];
    for (int row = 0; row < 4; row++) {
        tmp[row] = m->m[0*4 + row] * v[0]
                 + m->m[1*4 + row] * v[1]
                 + m->m[2*4 + row] * v[2]
                 + m->m[3*4 + row] * v[3];
    }
    out[0] = tmp[0]; out[1] = tmp[1]; out[2] = tmp[2]; out[3] = tmp[3];
}

static inline void m4_perspective(m4 *out, float fov_rad,
                                  float aspect, float near, float far) {
    float f = 1.0f / tanf(fov_rad * 0.5f);
    for (int i = 0; i < 16; i++) out->m[i] = 0.0f;
    out->m[0]  = f / aspect;
    out->m[5]  = f;
    out->m[10] = (far + near) / (near - far);
    out->m[11] = -1.0f;
    out->m[14] = (2.0f * far * near) / (near - far);
}

static inline void m4_perspective_vk(m4 *out, float fov_rad,
                                     float aspect, float near, float far,
                                     int flip_y) {
    float f = 1.0f / tanf(fov_rad * 0.5f);
    for (int i = 0; i < 16; i++) out->m[i] = 0.0f;
    out->m[0]  = f / aspect;
    out->m[5]  = flip_y ? -f : f;
    out->m[10] = far / (near - far);
    out->m[11] = -1.0f;
    out->m[14] = (far * near) / (near - far);
}

static inline void m4_look_at(m4 *out, vec3 eye, vec3 center, vec3 up) {
    vec3 f = vnormalize(vsub(center, eye));
    vec3 r = vnormalize(vcross(f, up));
    vec3 u = vcross(r, f);

    m4_identity(out);
    out->m[0] = r.x;  out->m[4] = r.y;  out->m[8]  = r.z;
    out->m[1] = u.x;  out->m[5] = u.y;  out->m[9]  = u.z;
    out->m[2] = -f.x; out->m[6] = -f.y; out->m[10] = -f.z;
    out->m[12] = -vdot(r, eye);
    out->m[13] = -vdot(u, eye);
    out->m[14] =  vdot(f, eye);
}

static inline void m4_ortho(m4 *out,
                            float left, float right,
                            float bottom, float top,
                            float near, float far) {
    m4_identity(out);
    out->m[0]  = 2.0f / (right - left);
    out->m[5]  = 2.0f / (top - bottom);
    out->m[10] = 2.0f / (near - far);
    out->m[12] = (right + left) / (right - left);
    out->m[13] = (top + bottom) / (top - bottom);
    out->m[14] = (near + far) / (near - far);
}

static inline void m4_ortho_vk(m4 *out,
                               float left, float right,
                               float bottom, float top,
                               float near, float far,
                               int flip_y) {
    m4_identity(out);
    out->m[0]  = 2.0f / (right - left);
    out->m[5]  = flip_y ? -2.0f / (top - bottom) : 2.0f / (top - bottom);
    out->m[10] = 1.0f / (near - far);
    out->m[12] = -(right + left) / (right - left);
    out->m[13] = flip_y ? -(top + bottom) / (top - bottom)
                        :  (top + bottom) / (top - bottom);
    out->m[14] = near / (near - far);
}

static inline void m4_transpose(m4 *out, const m4 *a) {
    m4 tmp;
    for (int col = 0; col < 4; col++)
        for (int row = 0; row < 4; row++)
            tmp.m[col * 4 + row] = a->m[row * 4 + col];
    *out = tmp;
}

static inline int m4_inverse(m4 *out, const m4 *a) {
    const float *m = a->m;
    float inv[16];

    inv[0] =  m[5]*m[10]*m[15] - m[5]*m[11]*m[14] - m[9]*m[6]*m[15]
            + m[9]*m[7]*m[14] + m[13]*m[6]*m[11] - m[13]*m[7]*m[10];
    inv[4] = -m[4]*m[10]*m[15] + m[4]*m[11]*m[14] + m[8]*m[6]*m[15]
            - m[8]*m[7]*m[14] - m[12]*m[6]*m[11] + m[12]*m[7]*m[10];
    inv[8] =  m[4]*m[9]*m[15] - m[4]*m[11]*m[13] - m[8]*m[5]*m[15]
            + m[8]*m[7]*m[13] + m[12]*m[5]*m[11] - m[12]*m[7]*m[9];
    inv[12] = -m[4]*m[9]*m[14] + m[4]*m[10]*m[13] + m[8]*m[5]*m[14]
            - m[8]*m[6]*m[13] - m[12]*m[5]*m[10] + m[12]*m[6]*m[9];
    inv[1] = -m[1]*m[10]*m[15] + m[1]*m[11]*m[14] + m[9]*m[2]*m[15]
            - m[9]*m[3]*m[14] - m[13]*m[2]*m[11] + m[13]*m[3]*m[10];
    inv[5] =  m[0]*m[10]*m[15] - m[0]*m[11]*m[14] - m[8]*m[2]*m[15]
            + m[8]*m[3]*m[14] + m[12]*m[2]*m[11] - m[12]*m[3]*m[10];
    inv[9] = -m[0]*m[9]*m[15] + m[0]*m[11]*m[13] + m[8]*m[1]*m[15]
            - m[8]*m[3]*m[13] - m[12]*m[1]*m[11] + m[12]*m[3]*m[9];
    inv[13] = m[0]*m[9]*m[14] - m[0]*m[10]*m[13] - m[8]*m[1]*m[14]
            + m[8]*m[2]*m[13] + m[12]*m[1]*m[10] - m[12]*m[2]*m[9];
    inv[2] =  m[1]*m[6]*m[15] - m[1]*m[7]*m[14] - m[5]*m[2]*m[15]
            + m[5]*m[3]*m[14] + m[13]*m[2]*m[7] - m[13]*m[3]*m[6];
    inv[6] = -m[0]*m[6]*m[15] + m[0]*m[7]*m[14] + m[4]*m[2]*m[15]
            - m[4]*m[3]*m[14] - m[12]*m[2]*m[7] + m[12]*m[3]*m[6];
    inv[10] = m[0]*m[5]*m[15] - m[0]*m[7]*m[13] - m[4]*m[1]*m[15]
            + m[4]*m[3]*m[13] + m[12]*m[1]*m[7] - m[12]*m[3]*m[5];
    inv[14] = -m[0]*m[5]*m[14] + m[0]*m[6]*m[13] + m[4]*m[1]*m[14]
            - m[4]*m[2]*m[13] - m[12]*m[1]*m[6] + m[12]*m[2]*m[5];
    inv[3] = -m[1]*m[6]*m[11] + m[1]*m[7]*m[10] + m[5]*m[2]*m[11]
            - m[5]*m[3]*m[10] - m[9]*m[2]*m[7] + m[9]*m[3]*m[6];
    inv[7] =  m[0]*m[6]*m[11] - m[0]*m[7]*m[10] - m[4]*m[2]*m[11]
            + m[4]*m[3]*m[10] + m[8]*m[2]*m[7] - m[8]*m[3]*m[6];
    inv[11] = -m[0]*m[5]*m[11] + m[0]*m[7]*m[9] + m[4]*m[1]*m[11]
            - m[4]*m[3]*m[9] - m[8]*m[1]*m[7] + m[8]*m[3]*m[5];
    inv[15] = m[0]*m[5]*m[10] - m[0]*m[6]*m[9] - m[4]*m[1]*m[10]
            + m[4]*m[2]*m[9] + m[8]*m[1]*m[6] - m[8]*m[2]*m[5];

    float det = m[0]*inv[0] + m[1]*inv[4] + m[2]*inv[8] + m[3]*inv[12];
    if (fabsf(det) < 1e-8f) return 0;

    float invDet = 1.0f / det;
    for (int i = 0; i < 16; i++) out->m[i] = inv[i] * invDet;
    return 1;
}

static inline vec3 m4_transform_point(const m4 *m, vec3 v) {
    float in[4] = { v.x, v.y, v.z, 1.0f };
    float out[4];
    m4_mul_vec4(out, m, in);
    if (fabsf(out[3]) < 1e-8f)
        return vmake(out[0], out[1], out[2]);
    float inv_w = 1.0f / out[3];
    return vmake(out[0] * inv_w, out[1] * inv_w, out[2] * inv_w);
}

static inline vec3 m4_transform_dir(const m4 *m, vec3 v) {
    float in[4] = { v.x, v.y, v.z, 0.0f };
    float out[4];
    m4_mul_vec4(out, m, in);
    return vmake(out[0], out[1], out[2]);
}

//--++=== quaternions(хз правильно ли написал? похуй) ===++--

static inline quat qidentity(void) {
    quat q = {0, 0, 0, 1};
    return q;
}

static inline quat qaxis(vec3 axis, float angle) {
    vec3 a = vnormalize(axis);
    float half = angle * 0.5f;
    float s = sinf(half);
    quat q = {a.x * s, a.y * s, a.z * s, cosf(half)};
    return q;
}

static inline quat qmul(quat a, quat b) {
    quat q = {
        a.w*b.x + a.x*b.w + a.y*b.z - a.z*b.y,
        a.w*b.y - a.x*b.z + a.y*b.w + a.z*b.x,
        a.w*b.z + a.x*b.y - a.y*b.x + a.z*b.w,
        a.w*b.w - a.x*b.x - a.y*b.y - a.z*b.z
    };
    return q;
}

static inline quat qnormalize(quat q) {
    float len = sqrtf(q.x*q.x + q.y*q.y + q.z*q.z + q.w*q.w);
    if (len < 1e-8f) return qidentity();
    float inv = 1.0f / len;
    quat r = {q.x*inv, q.y*inv, q.z*inv, q.w*inv};
    return r;
}

static inline void q_to_m4(m4 *out, quat q) {
    q = qnormalize(q);
    float xx = q.x*q.x, yy = q.y*q.y, zz = q.z*q.z;
    float xy = q.x*q.y, xz = q.x*q.z, yz = q.y*q.z;
    float wx = q.w*q.x, wy = q.w*q.y, wz = q.w*q.z;

    m4_identity(out);
    out->m[0]  = 1.0f - 2.0f*(yy + zz);
    out->m[1]  = 2.0f*(xy + wz);
    out->m[2]  = 2.0f*(xz - wy);
    out->m[4]  = 2.0f*(xy - wz);
    out->m[5]  = 1.0f - 2.0f*(xx + zz);
    out->m[6]  = 2.0f*(yz + wx);
    out->m[8]  = 2.0f*(xz + wy);
    out->m[9]  = 2.0f*(yz - wx);
    out->m[10] = 1.0f - 2.0f*(xx + yy);
}

static inline quat qslerp(quat a, quat b, float t) {
    float dot = a.x*b.x + a.y*b.y + a.z*b.z + a.w*b.w;
    if (dot < 0.0f) { b.x = -b.x; b.y = -b.y; b.z = -b.z; b.w = -b.w; dot = -dot; }
    if (dot > 1.0f) dot = 1.0f;
    if (dot < -1.0f) dot = -1.0f;
    if (dot > 0.9995f) {
        quat r = {
            a.x + t*(b.x-a.x), a.y + t*(b.y-a.y),
            a.z + t*(b.z-a.z), a.w + t*(b.w-a.w)
        };
        return qnormalize(r);
    }
    float theta = acosf(dot);
    float sinTheta = sinf(theta);
    float wa = sinf((1.0f - t) * theta) / sinTheta;
    float wb = sinf(t * theta) / sinTheta;
    quat r = {
        a.x*wa + b.x*wb, a.y*wa + b.y*wb,
        a.z*wa + b.z*wb, a.w*wa + b.w*wb
    };
    return r;
}

//--++=== matrix 3x3 ===++--

static inline void m3_identity(m3 *out) {
    for (int i = 0; i < 9; i++) out->m[i] = 0.0f;
    out->m[0] = out->m[4] = out->m[8] = 1.0f;
}

static inline void m3_from_m4(m3 *out, const m4 *m) {
    out->m[0] = m->m[0]; out->m[3] = m->m[4]; out->m[6] = m->m[8];
    out->m[1] = m->m[1]; out->m[4] = m->m[5]; out->m[7] = m->m[9];
    out->m[2] = m->m[2]; out->m[5] = m->m[6]; out->m[8] = m->m[10];
}

static inline void m3_transpose(m3 *out, const m3 *a) {
    m3 tmp;
    tmp.m[0] = a->m[0]; tmp.m[1] = a->m[3]; tmp.m[2] = a->m[6];
    tmp.m[3] = a->m[1]; tmp.m[4] = a->m[4]; tmp.m[5] = a->m[7];
    tmp.m[6] = a->m[2]; tmp.m[7] = a->m[5]; tmp.m[8] = a->m[8];
    *out = tmp;
}

static inline int m3_inverse(m3 *out, const m3 *a) {
    const float *m = a->m;
    float det = m[0] * (m[4]*m[8] - m[7]*m[5])
              - m[3] * (m[1]*m[8] - m[7]*m[2])
              + m[6] * (m[1]*m[5] - m[4]*m[2]);

    if (fabsf(det) < 1e-8f) return 0;

    float inv = 1.0f / det;
    m3 r;
    r.m[0] =  (m[4]*m[8] - m[7]*m[5]) * inv;
    r.m[1] = -(m[1]*m[8] - m[7]*m[2]) * inv;
    r.m[2] =  (m[1]*m[5] - m[4]*m[2]) * inv;
    r.m[3] = -(m[3]*m[8] - m[6]*m[5]) * inv;
    r.m[4] =  (m[0]*m[8] - m[6]*m[2]) * inv;
    r.m[5] = -(m[0]*m[5] - m[3]*m[2]) * inv;
    r.m[6] =  (m[3]*m[7] - m[6]*m[4]) * inv;
    r.m[7] = -(m[0]*m[7] - m[6]*m[1]) * inv;
    r.m[8] =  (m[0]*m[4] - m[3]*m[1]) * inv;
    *out = r;
    return 1;
}

static inline int m3_normal_matrix(m3 *out, const m4 *model) {
    m3 tmp, inv;
    m3_from_m4(&tmp, model);
    if (!m3_inverse(&inv, &tmp)) {
        m3_identity(out);
        return 0;
    }
    m3_transpose(out, &inv);
    return 1;
}

//--++=== vector (x, y) ===++--

static inline vec2 v2make(float x, float y) {
    vec2 v = {x, y};
    return v;
}

static inline vec2 v2add(vec2 a, vec2 b) {
    return v2make(a.x + b.x, a.y + b.y);
}

static inline vec2 v2sub(vec2 a, vec2 b) {
    return v2make(a.x - b.x, a.y - b.y);
}

static inline vec2 v2scale(vec2 v, float s) {
    return v2make(v.x * s, v.y * s);
}

static inline float v2dot(vec2 a, vec2 b) {
    return a.x * b.x + a.y * b.y;
}

static inline float v2cross(vec2 a, vec2 b) {
    return a.x * b.y - a.y * b.x;
}

static inline float v2length_sq(vec2 v) {
    return v2dot(v, v);
}

static inline float v2length(vec2 v) {
    return sqrtf(v2length_sq(v));
}

static inline vec2 v2normalize(vec2 v) {
    float len = v2length(v);
    if (len < 1e-8f) return v2make(0.0f, 0.0f);
    return v2scale(v, 1.0f / len);
}

static inline vec2 v2lerp(vec2 a, vec2 b, float t) {
    return v2make(a.x + t * (b.x - a.x),
                  a.y + t * (b.y - a.y));
}

static inline vec2 v2perp(vec2 v) {
    return v2make(-v.y, v.x);
}

static inline vec2 v2rotate(vec2 v, float angle) {
    float c = cosf(angle);
    float s = sinf(angle);
    return v2make(v.x * c - v.y * s,
                  v.x * s + v.y * c);
}

//--++=== vector (x, y, z, w) ===++--

static inline vec4 v4make(float x, float y, float z, float w) {
    vec4 v = {x, y, z, w};
    return v;
}

static inline vec4 v4add(vec4 a, vec4 b) {
    return v4make(a.x + b.x, a.y + b.y, a.z + b.z, a.w + b.w);
}

static inline vec4 v4sub(vec4 a, vec4 b) {
    return v4make(a.x - b.x, a.y - b.y, a.z - b.z, a.w - b.w);
}

static inline vec4 v4scale(vec4 v, float s) {
    return v4make(v.x * s, v.y * s, v.z * s, v.w * s);
}

static inline vec4 v4mul(vec4 a, vec4 b) {
    return v4make(a.x * b.x, a.y * b.y, a.z * b.z, a.w * b.w);
}

static inline float v4dot(vec4 a, vec4 b) {
    return a.x * b.x + a.y * b.y + a.z * b.z + a.w * b.w;
}

static inline float v4length_sq(vec4 v) {
    return v4dot(v, v);
}

static inline float v4length(vec4 v) {
    return sqrtf(v4length_sq(v));
}

static inline vec4 v4normalize(vec4 v) {
    float len = v4length(v);
    if (len < 1e-8f) return v4make(0.0f, 0.0f, 0.0f, 0.0f);
    return v4scale(v, 1.0f / len);
}

static inline vec4 v4lerp(vec4 a, vec4 b, float t) {
    return v4make(a.x + t * (b.x - a.x),
                  a.y + t * (b.y - a.y),
                  a.z + t * (b.z - a.z),
                  a.w + t * (b.w - a.w));
}

static inline vec4 v4clamp(vec4 v, float lo, float hi) {
    return v4make(
        v.x < lo ? lo : (v.x > hi ? hi : v.x),
        v.y < lo ? lo : (v.y > hi ? hi : v.y),
        v.z < lo ? lo : (v.z > hi ? hi : v.z),
        v.w < lo ? lo : (v.w > hi ? hi : v.w)
    );
}

static inline vec4 v4divide_by_w(vec4 v) {
    if (fabsf(v.w) < 1e-8f)
        return v4make(v.x, v.y, v.z, 0.0f);
    float inv = 1.0f / v.w;
    return v4make(v.x * inv, v.y * inv, v.z * inv, 1.0f);
}

static inline vec4 v4from_vec3(vec3 v, float w) {
    return v4make(v.x, v.y, v.z, w);
}

static inline vec3 v4to_vec3(vec4 v) {
    return vmake(v.x, v.y, v.z);
}

//--++=== matrix 2x2 ===++--

static inline void m2_identity(m2 *out) {
    out->m[0] = 1.0f; out->m[1] = 0.0f;
    out->m[2] = 0.0f; out->m[3] = 1.0f;
}

static inline void m2_multiply(m2 *out, const m2 *a, const m2 *b) {
    m2 tmp;
    for (int col = 0; col < 2; col++) {
        for (int row = 0; row < 2; row++) {
            float sum = 0.0f;
            for (int k = 0; k < 2; k++)
                sum += a->m[k * 2 + row] * b->m[col * 2 + k];
            tmp.m[col * 2 + row] = sum;
        }
    }
    *out = tmp;
}

static inline void m2_transpose(m2 *out, const m2 *a) {
    m2 tmp;
    tmp.m[0] = a->m[0]; tmp.m[1] = a->m[2];
    tmp.m[2] = a->m[1]; tmp.m[3] = a->m[3];
    *out = tmp;
}

static inline float m2_determinant(const m2 *a) {
    return a->m[0] * a->m[3] - a->m[2] * a->m[1];
}

static inline int m2_inverse(m2 *out, const m2 *a) {
    float det = m2_determinant(a);
    if (fabsf(det) < 1e-8f) return 0;

    float inv_det = 1.0f / det;
    m2 tmp;
    tmp.m[0] =  a->m[3] * inv_det;
    tmp.m[1] = -a->m[1] * inv_det;
    tmp.m[2] = -a->m[2] * inv_det;
    tmp.m[3] =  a->m[0] * inv_det;

    *out = tmp;
    return 1;
}

static inline void m2_scale(m2 *out, float sx, float sy) {
    m2_identity(out);
    out->m[0] = sx;
    out->m[3] = sy;
}

static inline void m2_rotate(m2 *out, float angle) {
    float c = cosf(angle);
    float s = sinf(angle);
    out->m[0] =  c; out->m[1] =  s;
    out->m[2] = -s; out->m[3] =  c;
}

static inline void m2_from_angle(m2 *out, float angle) {
    m2_rotate(out, angle);
}

static inline vec2 m2_mul_vec2(const m2 *m, vec2 v) {
    return v2make(
        m->m[0] * v.x + m->m[2] * v.y,
        m->m[1] * v.x + m->m[3] * v.y
    );
}

//--++=== quat/vec helpers ===++--

static inline vec3 q_rotate_vec3(quat q, vec3 v) {
    vec3 qv = vmake(q.x, q.y, q.z);
    vec3 t = vscale(vcross(qv, v), 2.0f);
    return vadd(v, vadd(vscale(t, q.w), vcross(qv, t)));
}

static inline vec3 vrotate_axis(vec3 v, vec3 axis, float angle) {
    quat q = qaxis(axis, angle);
    return q_rotate_vec3(q, v);
}

static inline void m4_rotate_quat(m4 *out, quat q) {
    q_to_m4(out, q);
}

static inline void m4_rotate_axis(m4 *out, vec3 axis, float angle) {
    quat q = qaxis(axis, angle);
    q_to_m4(out, q);
}

#endif
