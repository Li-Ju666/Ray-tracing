#pragma once

#define GLM_FORCE_RADIANS
#include <glm/glm.hpp>
#include <glm/gtc/matrix_transform.hpp>
#include <glm/gtx/component_wise.hpp>
#include <iostream>
#include <math.h>

namespace rt {

glm::vec3 random_in_unit_sphere(); 
glm::vec3 reflect(const glm::vec3&, const glm::vec3&); 
bool refract(const glm::vec3&, const glm::vec3&, float, glm::vec3&); 
float schlick(float, float); 

class material; 
class HitRecord; 
class lambertian; 
class metal; 

class Hitable {
public:
    virtual bool hit(const Ray& r, float t_min, float t_max, HitRecord& rec) const = 0;
    material* mat_ptr;
};


struct HitRecord {
    float t;
    glm::vec3 p;
    glm::vec3 normal;
    material* mat_ptr; 
};

class material {
    public: 
        virtual bool scatter(const Ray& r_in, const HitRecord& rec, 
            glm::vec3& attenuation, Ray& scattered) const = 0; 
        virtual glm::vec3 emit(const HitRecord& rec) {
            return glm::vec3(0.0); 
        }
    public:
        glm::vec3 albedo; 
};

class lambertian : public material {
    public: 
        lambertian(const glm::vec3& a) {albedo = a; 
        }
        virtual bool scatter(const Ray& r_in, const HitRecord& rec, 
            glm::vec3& attenuation, Ray& scattered) const {
            //std::cout << "I am scattering! " << std::endl; 
            //_sleep(1000); 
            glm::vec3 target = rec.p + rec.normal + random_in_unit_sphere(); 
            scattered = Ray(rec.p, target - rec.p); 
            attenuation = albedo; 
            return true; 
        }
};

class metal : public material {
    public: 
        metal(const glm::vec3& a) { albedo = a; fuzz = 0; }
        metal(const glm::vec3& a, float f) { albedo = a; fuzz = f; }

        virtual bool scatter(const Ray& r_in, const HitRecord& rec,
            glm::vec3& attenuation, Ray& scattered) const {
            glm::vec3 reflected = reflect(glm::normalize(r_in.direction()), rec.normal); 
            scattered = Ray(rec.p, reflected+fuzz*random_in_unit_sphere()); 
            attenuation = albedo; 
            return (glm::dot(scattered.direction(), rec.normal) > 0); 
        }
        float fuzz; 
};

class dielectric : public material {
public: 
    dielectric(float ri) { ref_idx = ri;  }
    virtual bool scatter(const Ray& r_in, const HitRecord& rec, glm::vec3& attenuation, Ray& scattered) const {
        glm::vec3 outward_normal; 
        glm::vec3 reflected = reflect(r_in.direction(), rec.normal); 
        float ni_over_nt; 
        attenuation = glm::vec3(1.0f, 1.0f, 1.0f); 
        glm::vec3 refracted; 
        float reflect_prob; 
        float cosine; 
        if (glm::dot(r_in.direction(), rec.normal) > 0) {
            outward_normal = -rec.normal; 
            ni_over_nt = ref_idx; 
            cosine = ref_idx * glm::dot(r_in.direction(), rec.normal) /
                glm::length(r_in.direction()); 
        }
        else {
            outward_normal = rec.normal; 
            ni_over_nt = 1.0 / ref_idx; 
            cosine = -glm::dot(r_in.direction(), rec.normal) /
                glm::length(r_in.direction()); 
        }
        if (refract(r_in.direction(), outward_normal, ni_over_nt, refracted)) {
            reflect_prob = schlick(cosine, ref_idx); 
            //scattered = Ray(rec.p, refracted); 
        }
        else {
            //scattered = Ray(rec.p, reflected); 
            reflect_prob = 1.0; 
            //return false; 
        }
        if ((float)rand() / RAND_MAX < reflect_prob) {
            scattered = Ray(rec.p, reflected); 
            //std::cout << "A " << std::endl; 
        }
        else {
            scattered = Ray(rec.p, refracted);
            //std::cout << "B " << std::endl;
        }
        return true; 

    }
    float ref_idx; 

};

class diffuse_light : public material {
public:
    diffuse_light(glm::vec3& light_color, float t) { 
        emit_light = light_color; 
        intensity = t; 
        //std::cout << "Emit light is " << emit[0] << std::endl; 
    }
    virtual bool scatter(const Ray& r_in, const HitRecord& rec, glm::vec3& attenuation, Ray& scattered) const {
        return false; 
    }

    virtual glm::vec3 emit(const HitRecord& rec) {
        return intensity*emit_light;
    }
    glm::vec3 emit_light;
    float intensity; 

};

} // namespace rt
