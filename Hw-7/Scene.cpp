//
// Created by Göksu Güvendiren on 2019-05-14.
//

#include "Scene.hpp"
#include <assert.h>

void Scene::buildBVH() {
    printf(" - Generating BVH...\n\n");
    this->bvh = new BVHAccel(objects, 1, BVHAccel::SplitMethod::NAIVE);
}

Intersection Scene::intersect(const Ray &ray) const
{
    return this->bvh->Intersect(ray);
}

void Scene::sampleLight(Intersection &pos, float &pdf) const//均匀采样场景中的一个光源
{
    float emit_area_sum = 0;
    for (uint32_t k = 0; k < objects.size(); ++k) {
        if (objects[k]->hasEmit()){
            emit_area_sum += objects[k]->getArea();
        }
    }
    float p = get_random_float() * emit_area_sum;
    emit_area_sum = 0;
    for (uint32_t k = 0; k < objects.size(); ++k) {
        if (objects[k]->hasEmit()){
            emit_area_sum += objects[k]->getArea();
            if (p <= emit_area_sum){
                objects[k]->Sample(pos, pdf);
                break;
            }
        }
    }
}

bool Scene::trace(
        const Ray &ray,
        const std::vector<Object*> &objects,
        float &tNear, uint32_t &index, Object **hitObject)
{
    *hitObject = nullptr;
    for (uint32_t k = 0; k < objects.size(); ++k) {
        float tNearK = kInfinity;
        uint32_t indexK;
        Vector2f uvK;
        if (objects[k]->intersect(ray, tNearK, indexK) && tNearK < tNear) {
            *hitObject = objects[k];
            tNear = tNearK;
            index = indexK;
        }
    }


    return (*hitObject != nullptr);
}

// Implementation of Path Tracing
const float EPS = 1e-3;
Vector3f Scene::castRay(const Ray &ray, int depth) const
{
    // TO DO Implement Path Tracing Algorithm here
    //首先计算光源的直接贡献
    Intersection direct_il;
    Intersection inter = intersect(ray);
    if(!inter.happened||inter.m==nullptr)return Vector3f(0.f);
    if (inter.m->hasEmission())
    {
        if (!depth)return inter.m->getEmission();//如果是直接看到光源，那就直接返回辐射度
        else return Vector3f(0.f);//如果不是直接看到光源，光源的贡献在上一层已经被计算过了
    }
    float PDF;
    sampleLight(direct_il, PDF);//采样了一个小面元
    Material* M = inter.m;
    const Vector3f& wi = ray.direction;
    const Vector3f& P = inter.coords;
    const Vector3f& N = inter.normal;
    Vector3f irradiance;
    Vector3f BRDF;
    //测试小面元光源能否直接照到该点
    Vector3f PtoL = direct_il.coords - P;
    float dist = PtoL.norm();
    if (dist != 0.0)
    {
        PtoL.x /= dist;
        PtoL.y /= dist;
        PtoL.z /= dist;
    }
    Ray test_ray(P, PtoL);
    Intersection test = intersect(test_ray);
    if (std::abs(test.distance - dist) > EPS)//不能直接照到
        irradiance = Vector3f(0.f);
    else
    {
        BRDF= M->eval(wi, PtoL, N);
        const Vector3f& NL = direct_il.normal;
        irradiance = direct_il.emit * BRDF * dotProduct(NL,-PtoL) * dotProduct(N, PtoL) / dist / dist / PDF;//计算光源在该点的Irradiance
    }
    const float P_RR = get_random_float();
    if (P_RR > RussianRoulette)return irradiance;
    Vector3f dir_out = M->sample(wi, N);
    Ray new_ray(P, dir_out);
    BRDF = M->eval(wi, dir_out, N);
    PDF = M->pdf(wi, dir_out, N);
    irradiance += castRay(new_ray, depth + 1) * BRDF * dotProduct(N, dir_out) / PDF / RussianRoulette;//Indirect Irradiance
    const float& inter_dist = inter.distance;
    return irradiance;
}