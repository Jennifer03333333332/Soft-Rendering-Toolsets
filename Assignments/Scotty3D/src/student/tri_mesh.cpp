
#include "../rays/tri_mesh.h"
#include "../rays/samplers.h"

namespace PT {

BBox Triangle::bbox() const {

    // TODO (PathTracer): Task 2 or 3
    // Compute the bounding box of the triangle.

    Vec3 p0 = vertex_list[v0].position;
    Vec3 p1 = vertex_list[v1].position;
    Vec3 p2 = vertex_list[v2].position;

    float minX = std::min({p0.x, p1.x, p2.x});
    float maxX = std::max({p0.x, p1.x, p2.x});
    float minY = std::min({p0.y, p1.y, p2.y});
    float maxY = std::max({p0.y, p1.y, p2.y});
    float minZ = std::min({p0.z, p1.z, p2.z});
    float maxZ = std::max({p0.z, p1.z, p2.z});
    // Beware of flat/zero-volume boxes! You may need to
    // account for that here, or later on in BBox::intersect.
    BBox box(Vec3(minX, minY, minZ), Vec3(maxX, maxY, maxZ)); // min max
    return box;
}

Trace Triangle::hit(const Ray& ray) const {

    // Each vertex contains a postion and surface normal
    Tri_Mesh_Vert v_0 = vertex_list[v0];
    Tri_Mesh_Vert v_1 = vertex_list[v1];
    Tri_Mesh_Vert v_2 = vertex_list[v2];
    // (void)v_0;
    // (void)v_1;
    // (void)v_2;

    // TODO (PathTracer): Task 2
    // Intersect the ray with the triangle defined by the three vertices.

    // Moller-Trumbore algorithm
    // if ray intersect with this triangle
    bool hit = true; 
    //Store u,v,t
    Vec3 uvt;
    float distance;

    Vec3 p0p1 = v_1.position - v_0.position;      // e1
    Vec3 p0p2 = v_2.position - v_0.position;      // e2
    Vec3 p0_rayOrigin = ray.point - v_0.position; // s
    float det = dot(cross(p0p1, ray.dir), p0p2) ;
    if(dot(cross(p0p1, ray.dir), p0p2) != 0) { // if denominator is zero
        Vec3 Temp_vec3(-1.0f * dot(cross(p0_rayOrigin, p0p2), ray.dir),
                       dot(cross(p0p1, ray.dir), p0_rayOrigin),
                       -1.0f * dot(cross(p0_rayOrigin, p0p2), p0p1));
        uvt = Temp_vec3 / (dot(cross(p0p1, ray.dir), p0p2));
        // How to judge if ray hits the triangle?
        
        // The hit point not in the triangle; t must >= 0
        if(uvt.x < 0 || uvt.y < 0 || (1.0f - uvt.x - uvt.y) < 0 || uvt.z < 0) {
            hit = false; 
        }
        // Also point should be within the ray.dist_bounds
        distance = std::abs((uvt.z * ray.dir).norm());
        if(distance < ray.dist_bounds.x || distance > ray.dist_bounds.y) {
            hit = false;
        }
    } else {
        // TODO
        // dot(cross(p0p1, ray.dir), p0p2) == 0
        // means cross(p0p1, ray.dir) == 0 or p0p2 == 0 or cross(p0p1, ray.dir)  vertical to p0p2 
        // (1) p0p1 X ray.dir == zero: ray is parallel to p0p1
        // (2) cross(p0p1, ray.dir) is vertical to p0p2: ray.dir is in this triangle's plane
        
        hit = false;//temporary

        // Vec3 normal_thistriangle = cross(p0p1,p0p2);
        // Vec3 p0_rayO = ray.point - v_0.position;
        // if(dot(p0_rayO,normal_thistriangle) == 0){//ray.point is in this plane

        // }
    }

    Trace ret;
    ret.origin = ray.point;
    ret.hit = hit;       // was there an intersection?
    if(hit){
        ret.distance = distance;   // at what distance did the intersection occur?
        ret.position = ray.at(uvt.z); // where was the intersection? o+td
        // what was the surface normal at the intersection?// (this should be interpolated between the three vertex normals)
        ret.normal = uvt.x*v_0.normal +
            uvt.y*v_1.normal +
            (1.0f - uvt.x - uvt.y)*v_2.normal;
    }

       
    return ret;
}

Triangle::Triangle(Tri_Mesh_Vert* verts, unsigned int v0, unsigned int v1, unsigned int v2)
    : vertex_list(verts), v0(v0), v1(v1), v2(v2) {
}

Vec3 Triangle::sample(Vec3 from) const {
    Tri_Mesh_Vert v_0 = vertex_list[v0];
    Tri_Mesh_Vert v_1 = vertex_list[v1];
    Tri_Mesh_Vert v_2 = vertex_list[v2];
    Samplers::Triangle sampler(v_0.position, v_1.position, v_2.position);
    Vec3 pos = sampler.sample();
    return (pos - from).unit();
}

float Triangle::pdf(Ray wray, const Mat4& T, const Mat4& iT) const {

    Ray tray = wray;
    tray.transform(iT);

    Trace trace = hit(tray);
    if(trace.hit) {
        trace.transform(T, iT.T());
        Vec3 v_0 = T * vertex_list[v0].position;
        Vec3 v_1 = T * vertex_list[v1].position;
        Vec3 v_2 = T * vertex_list[v2].position;
        float a = 2.0f / cross(v_1 - v_0, v_2 - v_0).norm();
        float g =
            (trace.position - wray.point).norm_squared() / std::abs(dot(trace.normal, wray.dir));
        return a * g;
    }
    return 0.0f;
}

void Tri_Mesh::build(const GL::Mesh& mesh, bool bvh) {

    use_bvh = bvh;
    verts.clear();
    triangle_bvh.clear();
    triangle_list.clear();

    for(const auto& v : mesh.verts()) {
        verts.push_back({v.pos, v.norm});
    }

    const auto& idxs = mesh.indices();

    std::vector<Triangle> tris;
    for(size_t i = 0; i < idxs.size(); i += 3) {
        tris.push_back(Triangle(verts.data(), idxs[i], idxs[i + 1], idxs[i + 2]));
    }

    if(use_bvh) {
        triangle_bvh.build(std::move(tris), 4);
    } else {
        triangle_list = List<Triangle>(std::move(tris));
    }
}

Tri_Mesh::Tri_Mesh(const GL::Mesh& mesh, bool use_bvh) {
    build(mesh, use_bvh);
}

Tri_Mesh Tri_Mesh::copy() const {
    Tri_Mesh ret;
    ret.verts = verts;
    ret.triangle_bvh = triangle_bvh.copy();
    ret.triangle_list = triangle_list.copy();
    ret.use_bvh = use_bvh;
    return ret;
}

BBox Tri_Mesh::bbox() const {
    if(use_bvh) return triangle_bvh.bbox();
    return triangle_list.bbox();
}

Trace Tri_Mesh::hit(const Ray& ray) const {
    if(use_bvh) return triangle_bvh.hit(ray);
    return triangle_list.hit(ray);
}

size_t Tri_Mesh::visualize(GL::Lines& lines, GL::Lines& active, size_t level,
                           const Mat4& trans) const {
    if(use_bvh) return triangle_bvh.visualize(lines, active, level, trans);
    return 0;
}

Vec3 Tri_Mesh::sample(Vec3 from) const {
    if(use_bvh) {
        die("Sampling BVH-based triangle meshes is not yet supported.");
    }
    return triangle_list.sample(from);
}

float Tri_Mesh::pdf(Ray ray, const Mat4& T, const Mat4& iT) const {
    if(use_bvh) {
        die("Sampling BVH-based triangle meshes is not yet supported.");
    }
    return triangle_list.pdf(ray, T, iT);
}

} // namespace PT
