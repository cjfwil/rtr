#include <stdio.h>


#pragma warning(push, 0)
#define STB_IMAGE_IMPLEMENTATION
#include "stb_image.h"

#define STB_IMAGE_WRITE_IMPLEMENTATION
#include "stb_image_write.h"
#pragma warning(pop)

typedef int unsigned u32;
typedef u32 b32;
typedef char unsigned u8;

#define Min(a, b) (a < b) ? a : b
#define Max(a, b) (a > b) ? a : b

#define Clamp(v, a, b) ((v < a) ? a : ((v > b) ? b : v))

#define true 1
#define false 0

#define PI 3.141592653589793238462f

struct v3
{
    float x;
    float y;
    float z;
};

struct matrix3x3 {
    float v[3*3];
};

v3 Vec3(float x, float y, float z)
{
    v3 result = {};
    result.x = x;
    result.y = y;
    result.z = z;
    return(result);
}


v3 operator/(v3 a, float b)
{
    v3 result = {};
    result.x = a.x / b;
    result.y = a.y / b;
    result.z = a.z / b;
    return(result);
}

v3 operator/(float a, v3 b)
{
    v3 result = {};
    result.x = a / b.x;
    result.y = a / b.y;
    result.z = a / b.z;
    return(result);
}

v3 operator*(v3 a, v3 b)
{
    v3 result = {};
    result.x = a.x * b.x;
    result.y = a.y * b.y;
    result.z = a.z * b.z;
    return(result);
}

v3 operator*(v3 a, float b)
{
    v3 result = {};
    result.x = a.x * b;
    result.y = a.y * b;
    result.z = a.z * b;
    return(result);
}

v3 operator*(float b, v3 a)
{
    v3 result = {};
    result = a * b;
    return(result);
}

v3 operator+(v3 a, v3 b)
{
    v3 result = {};
    result.x = a.x + b.x;
    result.y = a.y + b.y;
    result.z = a.z + b.z;
    return(result);
}

v3 operator-(v3 a, v3 b)
{
    v3 result = {};
    result.x = a.x - b.x;
    result.y = a.y - b.y;
    result.z = a.z - b.z;
    return(result);
}

float Dot(v3 a, v3 b)
{
    float result = 0.0f;
    result = a.x*b.x + a.y*b.y + a.z*b.z; 
    return(result);
}

v3 Cross(v3 a, v3 b)
{
    v3 v = {};
    v.x = a.y*b.z - a.z*b.y;
    v.y = a.z*b.x - a.x*b.z;
    v.z = a.x*b.y - a.y*b.x;
    return(v);
}

float Magnitude(v3 a)
{
    float result = 0.0f;
    result = sqrtf(Dot(a, a));
    return(result);
}

v3 Normalise(v3 a)
{
    v3 result = {};
    result = a / Magnitude(a);
    return(result);
}

v3 operator*(matrix3x3 m, v3 b)
{
    v3 result = {};
    result.x = Dot(Vec3(m.v[0], m.v[1], m.v[2]), b);
    result.y = Dot(Vec3(m.v[3], m.v[4], m.v[5]), b);
    result.z = Dot(Vec3(m.v[6], m.v[7], m.v[8]), b);
    return(result);
}

v3 operator*(v3 b, matrix3x3 m)
{
    v3 result = m * b;
    return(result);
}

matrix3x3 RotationMatrixX(float theta)
{
    matrix3x3 m = {
        1, 0, 0,
        0, cosf(theta), -sinf(theta),
        0, sinf(theta), cosf(theta),
    };
    return(m);
}

matrix3x3 RotationMatrixY(float theta)
{
    matrix3x3 m = {
        cosf(theta), 0, sinf(theta),
        0, 1, 0,
        -sinf(theta), 0, cosf(theta),
    };
    return(m);
}

matrix3x3 RotationMatrixZ(float theta)
{
    matrix3x3 m = {
        cosf(theta), -sinf(theta), 0,
        sinf(theta), cosf(theta), 0,
        0, 0, 1,
    };
    return(m);
}

union rgb {
    unsigned short data;
    struct {
        u8 r;
        u8 g;
        u8 b;
    };
};

struct primitive_plane {
    v3 normal;
    v3 origin;
};

struct primitive_sphere {
    v3 pos;
    float radius;
};

struct light {
    v3 pos;
};

float LinePlaneIntersection(v3 rayOrigin, v3 rayDir, primitive_plane p)
{
    float result = Dot(p.normal, (p.origin-rayOrigin)) / (Dot(p.normal, rayDir));
    return(result);
}

float LineSphereIntersection(v3 rayOrigin, v3 rayDir, primitive_sphere s)
{
    float result = INFINITY;
    float a = Dot(rayDir, rayDir);
    float b = (Dot(rayDir, rayOrigin-s.pos));
    float c = Dot(rayOrigin-s.pos, rayOrigin-s.pos) - (s.radius*s.radius);
    
    float discriminant = (b*b) - (a*c);
    if (discriminant > 0) {
        float t1 = (-b + sqrtf(discriminant)) / (a);
        float t2 = (-b - sqrtf(discriminant)) / (a);
        //TODO: account for inside of sphere
        if (t1 >= 0 && t2 >= 0) {
            result = Min(t1, t2);
        }
    }
    return(result);
}

enum SURFACE_TYPE {
    TEXTURE,
    REFRACTING,
    REFLECTING,
};

#define PLANE_NUM 6
#define SPHERE_NUM 5
enum PRIMITIVE_TYPE {
    SPHERE,
    PLANE,
};

#define PRIMITIVE_TOTAL PLANE_NUM + SPHERE_NUM
struct scene_primitive {
    PRIMITIVE_TYPE type;
    SURFACE_TYPE surfaceType;
    u32 index;
};

int main(void)
{
    u32 w = 1920;
    u32 h = 1080;
    u32 comp = 3;
    void* data = malloc(w*h*comp);
    for (u32 i = 0; i < w*h*comp; ++i) { ((u8*)data)[i] = 0; }
    
    float frustumW = 20.0f;
    float frustumH = frustumW * (h / (float)w);
    float frustumZ = 5.0f;
    
    v3 cameraPos = Vec3(0, 1, -2);
    v3 cameraDirection = Normalise(Vec3(0, 0, -1));
    
#define LIGHT_NUM 1
    light lights[LIGHT_NUM] = {};
    lights[0].pos = Vec3(0, 2, 0);
    //lights[1].pos = Vec3(0, 0, 0);
    
    float radii = 0.5f;
    
    
    SURFACE_TYPE spheresSurface[SPHERE_NUM] = { 
        REFLECTING, TEXTURE, REFLECTING, TEXTURE, TEXTURE
    };
    SURFACE_TYPE planesSurface[PLANE_NUM] = { 
        TEXTURE, TEXTURE, TEXTURE, TEXTURE, TEXTURE, TEXTURE, 
    };
    
    primitive_sphere spheres[SPHERE_NUM] = {};
    spheres[0].pos = Vec3(0, 2, 2);
    spheres[0].radius = 1.0f;
    spheres[1].pos = Vec3(1, 0, 1);
    spheres[1].radius = radii;
    spheres[2].pos = Vec3(2, 2, 1);
    spheres[2].radius = radii;
    spheres[3].pos = Vec3(-1, 0, 1);
    spheres[3].radius = radii;
    spheres[4].pos = Vec3(-2, 0, 1);
    spheres[4].radius = radii;
    
    int texW;
    int texH;
    int n = 3;
    
    unsigned char *texData = stbi_load("earth.jpg", &texW, &texH, &n, 4);
    
    int planeTexW;
    int planeTexH;
    unsigned char *planeTex = stbi_load("planeTex2.jpg", &planeTexW, &planeTexH, &n, 4);
    
    primitive_plane planes[PLANE_NUM] = {};
    planes[0].normal = Vec3(0, 1, 0);
    planes[0].origin = Vec3(0, -1.5f, 0);
    
    planes[1].normal = Vec3(0, 0, -1);
    planes[1].origin = Vec3(0, 0, 4);
    
    planes[2].normal = Vec3(0, 0, 1);
    planes[2].origin = Vec3(0, 0, -4);
    
    planes[3].normal = Vec3(-1, 0, 0);
    planes[3].origin = Vec3(4, 0, 0);
    
    planes[4].normal = Vec3(1, 0, 0);
    planes[4].origin = Vec3(-4, 0, 0);
    
    planes[5].normal = Vec3(0, -1, 0);
    planes[5].origin = Vec3(0, 4, 0);
    
    
    scene_primitive scenePrimitives[PRIMITIVE_TOTAL] = {};
    
    u32 scenePrimitiveIndex = 0;
    for (u32 i = 0; i < SPHERE_NUM; ++i) {
        scene_primitive s = {};
        s.type = SPHERE;
        s.index = i;
        s.surfaceType = spheresSurface[i];
        scenePrimitives[scenePrimitiveIndex++] = s;
    }
    for (u32 i = 0; i < PLANE_NUM; ++i) {
        scene_primitive s = {};
        s.type = PLANE;
        s.index = i;
        s.surfaceType = planesSurface[i];
        scenePrimitives[scenePrimitiveIndex++] = s;
    }
    
    for (u32 x = 0; x < w; ++x) {
        for (u32 y = 0; y < h; ++y) {
            rgb finalClr = {};
            
            float px = ((x / (float)w) - 0.5f) * frustumW;
            float py = ((y / (float)h) - 0.5f) * frustumH;
            
            float nearestT = INFINITY;
            
            float brightness = 0.0f;
            float ambientBrightness = 0.3f;
            
            v3 rayOrigin = cameraPos;
            v3 rayDir = Normalise(Vec3(px, py, frustumZ) * RotationMatrixY(0));
            
            rgb pixelColour = {};
            pixelColour.r = 0xff;
            pixelColour.g = 0xff;
            pixelColour.b = 0xff;
            
            for (int i = 0; i < PRIMITIVE_TOTAL; ++i) {
                int currentPrimitiveIndex = i;
                v3 normal = {};
                v3 lightRayDir = {};
                v3 hitPoint = {};
                scene_primitive primitive = scenePrimitives[i];
                b32 hit = false;
                if (primitive.type == SPHERE) {
                    primitive_sphere currentSphere = spheres[primitive.index];
                    float smallestT = LineSphereIntersection(rayOrigin, rayDir, currentSphere);
                    if (smallestT < nearestT && smallestT >= 0) {
                        nearestT = smallestT;
                        hit = true;
                        hitPoint = rayOrigin + (nearestT * (rayDir));
                        normal = Normalise(hitPoint - currentSphere.pos);
                    }
                    if (hit && primitive.surfaceType == REFLECTING) {
                        v3 sphereNormal = Normalise(hitPoint - currentSphere.pos);
                        v3 reflectionDir = Normalise(hitPoint - (2*Dot(hitPoint, sphereNormal)*sphereNormal));
                        float nearestReflectedT = INFINITY;
                        for (int j = 0; j < PRIMITIVE_TOTAL; ++j) {
                            if (currentPrimitiveIndex==j) continue;
                            scene_primitive reflectedPrimitive = scenePrimitives[j];
                            if (reflectedPrimitive.type == SPHERE) {
                                primitive_sphere reflectedSphere = spheres[reflectedPrimitive.index];
                                float t = LineSphereIntersection(hitPoint, reflectionDir, reflectedSphere);
                                if (t >= 0 && t < nearestReflectedT) {
                                    nearestReflectedT = t;
                                    hitPoint = hitPoint + (nearestReflectedT*reflectionDir);
                                    primitive = reflectedPrimitive;
                                    currentPrimitiveIndex = j;
                                }
                            } else if (reflectedPrimitive.type == PLANE) {
                                primitive_plane reflectedPlane = planes[reflectedPrimitive.index];
                                float t = LinePlaneIntersection(hitPoint, reflectionDir, reflectedPlane);
                                if (t >= 0 && t < nearestReflectedT) {
                                    nearestReflectedT = t;
                                    hitPoint = hitPoint + (nearestReflectedT*reflectionDir);
                                    primitive = reflectedPrimitive;
                                    currentPrimitiveIndex = j;
                                }
                            }
                            
                        }
                    }
                }
                if (primitive.type == PLANE) {
                    primitive_plane currentPlane = planes[primitive.index];
                    float t = LinePlaneIntersection(rayOrigin, rayDir, currentPlane);
                    if (t < nearestT && t >= 0) {
                        nearestT = t;
                        hit = true;
                        hitPoint = (rayOrigin + nearestT*(rayDir));
                        normal = currentPlane.normal;
                    }
                }
                
                b32 inShadow = false;
                if (hit) {
                    for (int l = 0; l < LIGHT_NUM; ++l) {
                        light currentLight = lights[l];
                        lightRayDir = Normalise(currentLight.pos - hitPoint);
                        
                        for (int j = 0; j < PRIMITIVE_TOTAL; ++j) {
                            if (currentPrimitiveIndex==j) continue;
                            scene_primitive shadowPrimitive = scenePrimitives[j];
                            
                            float shadowT = INFINITY;
                            
                            float deltaHitLight = Magnitude(currentLight.pos - hitPoint);
                            if (shadowPrimitive.type == SPHERE) {
                                primitive_sphere shadowSphere = spheres[shadowPrimitive.index];
                                shadowT = LineSphereIntersection(hitPoint, lightRayDir, shadowSphere);
                                
                            }
                            if (shadowPrimitive.type == PLANE) {
                                primitive_plane shadowPlane = planes[shadowPrimitive.index];
                                shadowT = LinePlaneIntersection(hitPoint, lightRayDir, shadowPlane);
                            }
                            
                            if (!(shadowT == INFINITY || deltaHitLight < shadowT) && shadowT >= 0.0f) {
                                inShadow = true;
                                brightness = ambientBrightness;
                                break;
                            }
                        }
                    }
                }
                
                if (hit && primitive.surfaceType == TEXTURE) {
                    if (primitive.type == SPHERE) {
                        primitive_sphere currentSphere = spheres[primitive.index];
                        v3 vp = Normalise(hitPoint - currentSphere.pos);
                        v3 vn = Vec3(0, 1, 0);
                        v3 ve = Vec3(1, 0, 0);
                        float phi = acosf(-Dot(vn, vp));
                        float v = phi / PI;
                        float a = Clamp(Dot(vp, ve) / sinf(phi), -1, 1);
                        float theta = acosf(a) / (2*PI);
                        float u = 1-theta;
                        if (Dot(Cross(vn, ve), vp) > 0) {
                            u = theta;
                        }
                        int tx = (1-u) * (texW-1);
                        int ty = (1-v) * (texH-1);
                        pixelColour = ((rgb*)texData)[tx + (ty * texW)];
                    }
                    if (primitive.type == PLANE) {
                        primitive_plane currentPlane = planes[primitive.index];
                        v3 a = Cross(currentPlane.normal, Vec3(1, 0, 0));
                        v3 b = Cross(currentPlane.normal, Vec3(0, 1, 0));
                        v3 c = Cross(currentPlane.normal, Vec3(0, 0, 1));
                        v3 maxAB = Dot(a, a) > Dot(b, b) ? a : b;
                        v3 ud = Normalise(Dot(maxAB, maxAB) > Dot(c, c) ? maxAB : c);
                        v3 vd = Cross(currentPlane.normal, ud);
                        float u = Dot(ud, hitPoint);
                        float v = Dot(vd, hitPoint);
                        int tx = (int)((1-u) * (planeTexW-1)) % planeTexW;
                        if (tx < 0) {
                            tx = planeTexW - abs(tx);
                        }
                        int ty = (int)((1-v) * (planeTexW-1)) % planeTexH;
                        if (ty < 0) {
                            ty = planeTexH - abs(ty);
                        }
                        pixelColour = ((rgb*)planeTex)[tx + (ty * planeTexW)];
                    }
                }
                
                if (hit && !inShadow) {
                    brightness = Clamp(Dot(normal, lightRayDir), 0, 1);
                }
            }
            
            
            
            finalClr.r = pixelColour.r * brightness;
            finalClr.g = pixelColour.g * brightness;
            finalClr.b = pixelColour.b * brightness;
            
            u32 dataIndex = (x + ((h-y-1) * w))*comp;
            ((u8*)data)[dataIndex++] = finalClr.r;
            ((u8*)data)[dataIndex++] = finalClr.g;
            ((u8*)data)[dataIndex++] = finalClr.b;
        }
    }
    
    stbi_write_bmp("image.bmp", w, h, comp, data);
    return(0);
}