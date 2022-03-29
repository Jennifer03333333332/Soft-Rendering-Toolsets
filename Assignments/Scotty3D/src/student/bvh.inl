
#include "../rays/bvh.h"
#include "debug.h"
#include <stack>

namespace PT {

// template<typename Primitive>
// int partition(std::vector<Primitive>& nums, int left, int right, int pivot_index) {
//     int pivot = nums[pivot_index].bbox().center();

//     std::swap(nums[pivot_index], nums[right]);
//     int store_index = left;

//     for (int i = left; i <= right; i++) {
//         if (nums[i].bbox().center() < pivot) {
//             std::swap(nums[store_index], nums[i]);
//             store_index++;
//         }
//     }
//     std::swap(nums[store_index], nums[right]);
//     return store_index;
// }

// template<typename Primitive>
// int quickselect(std::vector<Primitive>& nums, int left, int right, int k_smallest) {
//     if (left == right) return nums[left];

//     srand(time(0));
//     int pivot_index = left + (rand() % (right - left));
//     pivot_index = partition(nums, left, right, pivot_index);

//     if (k_smallest == pivot_index) return nums[k_smallest];

//     else if (k_smallest < pivot_index) return quickselect(nums, left, pivot_index - 1,
//     k_smallest);

//     return quickselect(nums, pivot_index + 1, right, k_smallest);
// }

// template<typename Primitive>
// bool comp_center(Primitive& prim1, Primitive& prim2, int axis){
//     return (prim1.bbox().center()[axis] < prim2.bbox().center()[axis]);
// }

struct Bucket {
    BBox bbox;
    int prim_count = 0;
};

struct Partition {
    BBox b_left, b_right;
    int left_prim_count = 0, right_prim_count = 0;
};

template<typename Primitive>
void BVH<Primitive>::build(std::vector<Primitive>&& prims, size_t max_leaf_size) {

    // NOTE (PathTracer):
    // This BVH is parameterized on the type of the primitive it contains. This allows
    // us to build a BVH over any type that defines a certain interface. Specifically,
    // we use this to both build a BVH over triangles within each Tri_Mesh, and over
    // a variety of Objects (which might be Tri_Meshes, Spheres, etc.) in Pathtracer.

    // The Primitive interface must implement these two functions:
    //      BBox bbox() const;
    //      Trace hit(const Ray& ray) const;
    // Hence, you may call bbox() and hit() on any value of type Primitive.

    // Finally, also note that while a BVH is a tree structure, our BVH nodes don't
    // contain pointers to children, but rather indicies. This is because instead
    // of allocating each node individually, the BVH class contains a vector that
    // holds all of the nodes. Hence, to get the child of a node, you have to
    // look up the child index in this vector (e.g. nodes[node.l]). Similarly,
    // to create a new node, don't allocate one yourself - use BVH::new_node, which
    // returns the index of a newly added node.

    // Keep these
    nodes.clear();
    primitives = std::move(prims);

    // TODO (PathTracer): Task 3
    // Construct a BVH from the given vector of primitives and maximum leaf
    // size configuration. The starter code builds a BVH with a
    // single leaf node (which is also the root) that encloses all the primitives.

    // Replace these

    // What I have: primitives, max_leaf_size
    const int bucket_number = 10;

    BBox box; // bbox for all the primitives
    for(const Primitive& prim : primitives) {
        box.enclose(prim.bbox());
    }
    nodes.push_back(new_node(box, 0, primitives.size(), 0, 0));//box,start index, size, left
    // child, right child
    root_idx = 0;

    // Build the tree in level order
    while(true) {//++root_idx
        if(nodes[root_idx].size <= max_leaf_size) return;

        // For axis x,y,z
        float best_cost = FLT_MAX;
        Vec3 best_cost_ofxyz(FLT_MAX, FLT_MAX, FLT_MAX);
        std::unordered_map<float, Partition> best_partition; // self-define
        Partition best_partition_xyz;
        for(int axis = 0; i < 3; i++) {
            Partition best_partition_thisaxis;
            // Sort the primitives in ascending order, override the compare using lambda
            std::sort(primitives.begin(), primitives.end(),
                      [](const Primitive& prim1, const Primitive& prim2) {
                          return (prim1.bbox().center()[axis] < prim2.bbox().center()[axis]);
                      });

            // Initialize buckets

            float bucket_interval = (box.max[axis] - box.min[axis]) / (float)bucket_number;
            for(float line = box.min[axis] + bucket_interval; line < box.max[axis]; line += bucket_interval) { // the interval is 1
                // range: box.min[axis] -> line; line -> box.max[axis];
                // from axis = line to (line+1 or box.max[axis](when out of range))
                Partition p;
                for(int index = 0; index < primitives.size(); index++) { 
                    // This primitive is on the left
                    if(primitives[index].bbox().center()[axis] < line) {
                        p.b_left.enclose(primitives[index].bbox());
                        p.left_prim_count++;
                    } else { // on the right
                        p.b_right.enclose(primitives[index].bbox());
                        p.right_prim_count++;
                    }
                }
                // end of one partition
                // SAH Cost(for this axis)
                float cost = b_left.bbox.surface_area() / box.surface_area() * (float)b_left.size +
                             b_right.bbox.surface_area() / box.surface_area() * (float)b_right.size;
                if(cost < best_cost_ofxyz[axis]) {
                    best_cost_ofxyz[axis] = cost;
                    best_partition_thisaxis = p;
                }
            }

            best_partition.insert({best_cost_ofxyz[axis], best_partition_thisaxis});

            // std::partition(primitives.begin(),primitives.end(),[](){ return bbox().center() < })
        }
        best_cost = std::min(best_cost_ofxyz[0], best_cost_ofxyz[1], best_cost_ofxyz[2]);
        if(best_partition[best_cost]) {
            best_partition_xyz = best_partition[best_cost];
            // old node
            nodes[root_idx].l = (size_t)nodes.size();
            nodes[root_idx].r = (size_t)nodes.size() + 1;
            // left
            nodes.push_back(
                new_node(best_partition_xyz.b_left, nodes.size(), best_partition_xyz.left_prim_count, 0, 0)
            );

            // right
            nodes.push_back(
                new_node(best_partition_xyz.b_right, nodes.size(), best_partition_xyz.right_prim_count,0, 0);
            );
            
        }
        root_idx++;
    }

    // stop rule
    //  while( nodes[nodes.size()-1].size > max_leaf_size){

    // }
}

template<typename Primitive> Trace BVH<Primitive>::hit(const Ray& ray) const {

    // TODO (PathTracer): Task 3
    // Implement ray - BVH intersection test. A ray intersects
    // with a BVH aggregate if and only if it intersects a primitive in
    // the BVH that is not an aggregate.

    // The starter code simply iterates through all the primitives.
    // Again, remember you can use hit() on any Primitive value.

    Trace ret;
    for(const Primitive& prim : primitives) {
        Trace hit = prim.hit(ray);
        ret = Trace::min(ret, hit);
    }
    return ret;
}

template<typename Primitive>
BVH<Primitive>::BVH(std::vector<Primitive>&& prims, size_t max_leaf_size) {
    build(std::move(prims), max_leaf_size);
}

template<typename Primitive> BVH<Primitive> BVH<Primitive>::copy() const {
    BVH<Primitive> ret;
    ret.nodes = nodes;
    ret.primitives = primitives;
    ret.root_idx = root_idx;
    return ret;
}

template<typename Primitive> bool BVH<Primitive>::Node::is_leaf() const {

    // A node is a leaf if l == r, since all interior nodes must have distinct children
    return l == r;
}

template<typename Primitive>
size_t BVH<Primitive>::new_node(BBox box, size_t start, size_t size, size_t l, size_t r) {
    Node n;
    n.bbox = box;
    n.start = start;
    n.size = size;
    n.l = l;
    n.r = r;
    nodes.push_back(n);
    return nodes.size() - 1;
}

template<typename Primitive> BBox BVH<Primitive>::bbox() const {
    return nodes[root_idx].bbox;
}

template<typename Primitive> std::vector<Primitive> BVH<Primitive>::destructure() {
    nodes.clear();
    return std::move(primitives);
}

template<typename Primitive> void BVH<Primitive>::clear() {
    nodes.clear();
    primitives.clear();
}

template<typename Primitive>
size_t BVH<Primitive>::visualize(GL::Lines& lines, GL::Lines& active, size_t level,
                                 const Mat4& trans) const {

    std::stack<std::pair<size_t, size_t>> tstack;
    tstack.push({root_idx, 0});
    size_t max_level = 0;

    if(nodes.empty()) return max_level;

    while(!tstack.empty()) {

        auto [idx, lvl] = tstack.top();
        max_level = std::max(max_level, lvl);
        const Node& node = nodes[idx];
        tstack.pop();

        Vec3 color = lvl == level ? Vec3(1.0f, 0.0f, 0.0f) : Vec3(1.0f);
        GL::Lines& add = lvl == level ? active : lines;

        BBox box = node.bbox;
        box.transform(trans);
        Vec3 min = box.min, max = box.max;

        auto edge = [&](Vec3 a, Vec3 b) { add.add(a, b, color); };

        edge(min, Vec3{max.x, min.y, min.z});
        edge(min, Vec3{min.x, max.y, min.z});
        edge(min, Vec3{min.x, min.y, max.z});
        edge(max, Vec3{min.x, max.y, max.z});
        edge(max, Vec3{max.x, min.y, max.z});
        edge(max, Vec3{max.x, max.y, min.z});
        edge(Vec3{min.x, max.y, min.z}, Vec3{max.x, max.y, min.z});
        edge(Vec3{min.x, max.y, min.z}, Vec3{min.x, max.y, max.z});
        edge(Vec3{min.x, min.y, max.z}, Vec3{max.x, min.y, max.z});
        edge(Vec3{min.x, min.y, max.z}, Vec3{min.x, max.y, max.z});
        edge(Vec3{max.x, min.y, min.z}, Vec3{max.x, max.y, min.z});
        edge(Vec3{max.x, min.y, min.z}, Vec3{max.x, min.y, max.z});

        if(!node.is_leaf()) {
            tstack.push({node.l, lvl + 1});
            tstack.push({node.r, lvl + 1});
        } else {
            for(size_t i = node.start; i < node.start + node.size; i++) {
                size_t c = primitives[i].visualize(lines, active, level - lvl, trans);
                max_level = std::max(c + lvl, max_level);
            }
        }
    }
    return max_level;
}

} // namespace PT
