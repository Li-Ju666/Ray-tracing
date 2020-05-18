#include "raytracing.h"
#include "ray.h"
#include "hitable.h"
#include "sphere.h"
#include "triangle.h"
#include "box.h"
//#include <openmpi.h>

#include "utils2.h"  // Used for OBJ-mesh loading
#include <stdlib.h>  // Needed for drand48()
std::string modelDir(void); 

namespace rt {

// Store scene (world) in a global variable for convenience
struct Scene {
    Sphere ground;
    std::vector<Sphere> spheres;
    //std::vector<material> materials; 
    //std::vector<Hitable> list; 
    std::vector<Box> boxes;
    std::vector<Triangle> mesh;
    Box mesh_bbox;
} g_scene;

bool hit_world(const Ray &r, float t_min, float t_max, HitRecord &rec)
{
    HitRecord temp_rec;
    bool hit_anything = false;
    float closest_so_far = t_max;

    if (g_scene.ground.hit(r, t_min, closest_so_far, temp_rec)) {
        hit_anything = true;
        closest_so_far = temp_rec.t;
        rec = temp_rec;
    }
    for (int i = 0; i < g_scene.spheres.size(); ++i) {
        //std::cout << "My material is" << g_scene.spheres[i].mat_ptr->albedo[0] <<
        //    g_scene.spheres[i].mat_ptr->albedo[1] << g_scene.spheres[i].mat_ptr->albedo[2] << std::endl;
        if (g_scene.spheres[i].hit(r, t_min, closest_so_far, temp_rec)) {
            hit_anything = true;
            closest_so_far = temp_rec.t;
            rec = temp_rec;
        }
    }
    for (int i = 0; i < g_scene.boxes.size(); ++i) {
        if (g_scene.boxes[i].hit(r, t_min, closest_so_far, temp_rec)) {
            hit_anything = true;
            closest_so_far = temp_rec.t;
            rec = temp_rec;
        }
    }
    for (int i = 0; i < g_scene.mesh.size(); ++i) {
        if (g_scene.mesh[i].hit(r, t_min, closest_so_far, temp_rec)) {
            hit_anything = true;
            closest_so_far = temp_rec.t;
            rec = temp_rec;
        }
    }
    return hit_anything;
}

glm::vec3 reflect(const glm::vec3& v, const glm::vec3& n) {
    return v - 2.0f * glm::dot(v, n) * n; 
}

bool refract(const glm::vec3& v, const glm::vec3& n, float ni_over_nt, glm::vec3& refracted) {
    glm::vec3 uv = glm::normalize(v); 
    float dt = glm::dot(uv, n); 
    float discriminant = 1.0 - ni_over_nt * ni_over_nt * (1 - dt * dt); 
    if (discriminant > 0) {
        refracted = ni_over_nt * (uv - n * dt) - n * sqrt(discriminant);
        return true;
    }
    else
        return false; 
}

float schlick(float cosine, float ref_idx) {
    float r0 = (1 - ref_idx) / (1 + ref_idx); 
    r0 = r0 * r0; 
    return r0 + (1 - r0) * pow((1 - cosine), 5); 
}
// This function should be called recursively (inside the function) for
// bouncing rays when you compute the lighting for materials, like this
//
// if (hit_world(...)) {
//     ...
//     return color(rtx, r_bounce, max_bounces - 1);
// }
//
// See Chapter 7 in the "Ray Tracing in a Weekend" book
glm::vec3 random_in_unit_sphere() {
    glm::vec3 p; 

    
    do {
        float randa = 2 * (float)rand() / RAND_MAX - 1;
        float randb = 2 * (float)rand() / RAND_MAX - 1;
        float randc = 2 * (float)rand() / RAND_MAX - 1;
        p = 2.0f * glm::vec3(randa, randb, randc)-
            glm::vec3(1.0f,1.0f,1.0f); 
        
    } while (glm::length(p) >= 1.0); 

    //std::cout << p[0] << ", " << p[1] << ", " << p[2] << std::endl; 
    return p; 
}
glm::vec3 color(RTContext& rtx, const Ray& r, int max_bounces)
{
    //std::cout << "Going to color the pixel! " << std::endl; 
    if (max_bounces < 0) return glm::vec3(0.0f);

    HitRecord rec;
    if (hit_world(r, 0.001f, 9999.0f, rec)) {
        rec.normal = glm::normalize(rec.normal);  // Always normalise before use!
        if (rtx.show_normals) {
            //std::cout << "Normal vector visualization! " << std::endl; 
            return rec.normal * 0.5f + 0.5f;
        }
        else {
            //std::cout << "going to color! " << std::endl; 
            Ray scattered; 
            glm::vec3 attenuation; 
            glm::vec3 emitted = rec.mat_ptr->emit(rec); 
            //std::cout << "going to scatter! " << std::endl; 
            //_sleep(1000); 
            if (rec.mat_ptr->scatter(r, rec, attenuation, scattered)) {

                //std::cout << "Scatter finished! " << std::endl; 
                //_sleep(1000);
                return emitted + attenuation*color(rtx, scattered, max_bounces - 1); 
            }
            else {
                //std::cout << "Emitted is " << emitted[0] << std::endl; 
                return emitted; 
            }
            // Implement lighting for materials here
            // ...
            //std::cout << "To color! " << std::endl; 
            //glm::vec3 target = rec.p + rec.normal + random_in_unit_sphere();
            //std::cout << "Going to lightning for material! " << std::endl;
            //return 0.5f * color(rtx, Ray(rec.p, target - rec.p), max_bounces - 1);
        }
        //return glm::vec3(0.0f);
    }

    else {// If no hit, return sky color
        glm::vec3 unit_direction = glm::normalize(r.direction());
        float t = 0.5f * (unit_direction.y + 1.0f);
        //return (1.0f - t) * rtx.ground_color + t * rtx.sky_color; 
        return glm::vec3(0.0f); 
    }
}

// MODIFY THIS FUNCTION!
void setupScene(RTContext &rtx, const char *filename)
{
    
    material* pinklamb = new lambertian(glm::vec3(0.8, 0.3, 0.3)); 
    material* yellowmet = new metal(glm::vec3(0.8, 0.6, 0.2), 0.05); 
    material* bluemet = new metal(glm::vec3(0.2, 0.6, 0.8), 0.05);
    material* redlight = new diffuse_light(glm::vec3(0.8f, 0.8f, 0.0f), 50); 
 //   material* mat3 = new metal(glm::vec3(0.8, 0.8, 0.8), 0.2); 
    material* glass = new dielectric(1.55f); 
    material* matbackground = new lambertian(glm::vec3(0.8, 0.8, 0.8)); 
    g_scene.ground = Sphere(glm::vec3(0.0f, -1000.5f, 0.0f), 1000.0f, matbackground);
    g_scene.spheres = {
        Sphere(glm::vec3(0.0f, 0.0f, 0.0f), 0.5f, yellowmet),
        Sphere(glm::vec3(1.2f, 0.0f, 0.0f), 0.5f, pinklamb),
        Sphere(glm::vec3(-1.2f, 0.0f, 0.0f), 0.5f, glass),
        
        Sphere(glm::vec3(-1.0f, 5.0f, 0.0f), 0.1, redlight), 
        //Sphere(glm::vec3(0.0f, 1.0f, 0.0f), 0.3, glass),

        Sphere(glm::vec3(0.0f, -0.3f, 1.0f), 0.2f, glass),
        Sphere(glm::vec3(0.3f, -0.3f, 0.6f), 0.2f, bluemet),
        Sphere(glm::vec3(-0.8f, -0.4f, -0.7f), 0.1f, yellowmet),

        //Sphere(glm::vec3(1.0f, 0.0f, 1.0f), 0.1, mat2),
        //Sphere(glm::vec3(-1.0f, 0.0f, -1.0f), 0.1f, mat1),
        //Sphere(glm::vec3(-1.0f, 0.0f, 1.5f), 0.1f, mat3),
        //Sphere(glm::vec3(0.0f, 0.0f, -1.0f), 0.1f, mat3),

    };

    /*g_scene.boxes = {
        //Box(glm::vec3(0.0f, -0.25f, 0.0f), glm::vec3(0.25f)),
        //Box(glm::vec3(1.0f, -0.25f, 0.0f), glm::vec3(0.25f)),
        //Box(glm::vec3(0.0f, 0.0f, 0.0f), glm::vec3(10.0f)),
    };

    OBJMesh mesh;
    objMeshLoad(mesh, modelDir()+"armadillo_lowpoly.obj");
    g_scene.mesh.clear();
    for (int i = 0; i < mesh.indices.size(); i += 3) {
        int i0 = mesh.indices[i + 0];
        int i1 = mesh.indices[i + 1];
        int i2 = mesh.indices[i + 2];
        glm::vec3 v0 = mesh.vertices[i0] + glm::vec3(0.0f, 0.135f, 0.0f);
        glm::vec3 v1 = mesh.vertices[i1] + glm::vec3(0.0f, 0.135f, 0.0f);
        glm::vec3 v2 = mesh.vertices[i2] + glm::vec3(0.0f, 0.135f, 0.0f);
        g_scene.mesh.push_back(Triangle(v0, v1, v2, bluemet));
    }*/
}

// MODIFY THIS FUNCTION!
void updateLine(RTContext &rtx, int y)
{
    int nx = rtx.width;
    int ny = rtx.height;
    float aspect = float(nx) / float(ny);
    int ns = 2; 
    glm::vec3 lower_left_corner(-1.0f * aspect, -1.0f, -1.0f);
    glm::vec3 horizontal(2.0f * aspect, 0.0f, 0.0f);
    glm::vec3 vertical(0.0f, 2.0f, 0.0f);
    glm::vec3 origin(0.0f, 0.0f, 0.0f);
    glm::mat4 world_from_view = glm::inverse(rtx.view);

    // You can try to parallelise this loop by uncommenting this line:
    #pragma omp parallel for schedule(dynamic)
    for (int x = 0; x < nx; ++x) {
        glm::vec3 c = glm::vec3(0.0f); 
        if (rtx.anti_alias) {
            for (int s = 0; s < ns; s++) {
                float u = (float(x) + (float)rand() / RAND_MAX) / float(nx);
                float v = (float(y) + (float)rand() / RAND_MAX) / float(ny);
                Ray r(origin, lower_left_corner + u * horizontal + v * vertical);
                r.A = glm::vec3(world_from_view * glm::vec4(r.A, 1.0f));
                r.B = glm::vec3(world_from_view * glm::vec4(r.B, 0.0f));
                c += color(rtx, r, rtx.max_bounces);
            }
            c /= float(ns); 
        }
        else {
            float u = float(x)/ float(nx);
            float v = float(y) / float(ny);
            Ray r(origin, lower_left_corner + u * horizontal + v * vertical);
            r.A = glm::vec3(world_from_view * glm::vec4(r.A, 1.0f));
            r.B = glm::vec3(world_from_view * glm::vec4(r.B, 0.0f));
            c = color(rtx, r, rtx.max_bounces);
        }
        if (rtx.current_frame <= 0) {
                glm::vec4 old = rtx.image[y*nx + x];
                rtx.image[y*nx + x] = glm::clamp(old / glm::max(1.0f, old.a), 0.0f, 1.0f);
            }
        //c = c / float(ns); 
        rtx.image[y*nx + x] += glm::vec4(c, 1.0f);
    }
}

void updateImage(RTContext &rtx)
{
    if (rtx.freeze) return;  // Skip update
    rtx.image.resize(rtx.width * rtx.height);  // Just in case...

    updateLine(rtx, rtx.current_line % rtx.height);

    if (rtx.current_frame < rtx.max_frames) {
        rtx.current_line += 1;
        if (rtx.current_line >= rtx.height) {
            rtx.current_frame += 1;
            rtx.current_line = rtx.current_line % rtx.height;
        }
    }
}

void resetImage(RTContext &rtx)
{
    rtx.image.clear();
    rtx.image.resize(rtx.width * rtx.height);
    rtx.current_frame = 0;
    rtx.current_line = 0;
    rtx.freeze = false;
}

void resetAccumulation(RTContext &rtx)
{
    rtx.current_frame = -1;
}

} // namespace rt
