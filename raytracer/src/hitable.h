#pragma once

#define GLM_FORCE_RADIANS
#include <glm/glm.hpp>
#include <glm/gtc/matrix_transform.hpp>
#include <glm/gtx/component_wise.hpp>
#include <iostream>

namespace rt {

glm::vec3 random_in_unit_sphere(); 
glm::vec3 reflect(const glm::vec3&, const glm::vec3&); 

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
        //glm::vec3 albedo; 
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
        //glm::vec3 albedo; 
        float fuzz; 
};


} // namespace rt
