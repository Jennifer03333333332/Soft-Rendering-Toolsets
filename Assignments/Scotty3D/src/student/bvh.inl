
#include "../rays/bvh.h"
#include "debug.h"
#include <stack>

namespace PT {

// struct Bucket {
//     BBox bbox;
//     int prim_count = 0;
// };

struct Partition {
    
    int axis = 0;
    BBox b_left = {} , b_right = {};
    int left_prim_count = 0, right_prim_count = 0;

    // Partition(){
    //     axis = 0;
    //     b_left = BBox(); b_right= BBox();
    // }
};

template<typename Primitive>
bool x_cmp(const Primitive& prim1, const Primitive& prim2) {  
    return (prim1.bbox().center().x < prim2.bbox().center().x);
}
template<typename Primitive>
bool y_cmp(const Primitive& prim1, const Primitive& prim2) {  
    return (prim1.bbox().center().y < prim2.bbox().center().y);
}
template<typename Primitive>
bool z_cmp(const Primitive& prim1, const Primitive& prim2) {  
    return (prim1.bbox().center().z < prim2.bbox().center().z);
}

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
    root_idx = 0;
    // TODO (PathTracer): Task 3
    // Construct a BVH from the given vector of primitives and maximum leaf
    // size configuration. The starter code builds a BVH with a
    // single leaf node (which is also the root) that encloses all the primitives.

    // Replace these

    // What I have: primitives, max_leaf_size


    //start
    const int bucket_number = 10;

    BBox box; // bbox for all the primitives
    for(const Primitive& prim : primitives) {
        box.enclose(prim.bbox());
    }
    new_node(box, 0, primitives.size(), 0, 0);//box,start index, size, left
    // child, right child
    int cur_idx = 0;

    // Build the tree in level order
    while(true) {//++root_idx
        if(nodes[root_idx].size <= max_leaf_size) return;
        //primitive start index and end index
        int start_index = (int)nodes[root_idx].start;
        int end_index = (int)nodes[root_idx].start + (int)nodes[root_idx].size;
        // For axis x,y,z
        float best_cost = FLT_MAX;
        Vec3 best_cost_ofxyz(FLT_MAX, FLT_MAX, FLT_MAX);
        std::vector< Partition> best_partition; // self-define
        Partition best_partition_xyz;
        for(int axis = 0; axis < 3; axis++) {
            Partition best_partition_thisaxis;
            // Sort the primitives in ascending order, override the compare using lambda
            //based on this sub list
            

            // Initialize buckets

            float bucket_interval = (nodes[root_idx].bbox.max[axis] - nodes[root_idx].bbox.min[axis]) / (float)bucket_number;
            for(float line = nodes[root_idx].bbox.min[axis] + bucket_interval; line < nodes[root_idx].bbox.max[axis]; line += bucket_interval) { // the interval is 1
                // range: bbox.min[axis] -> line; line -> bbox.max[axis];
                // from axis = line to (line+1 or box.max[axis](when out of range))
                //it : the first one of the right part
                auto it = std::partition(primitives.begin(), primitives.end(), [](const Primitives&a){
                    return a.bbox().center() < line;
                });
                int index_intervel = it - primitives.begin();
                Partition p; p.axis = axis;
                for(int index = start_index; index < end_index; index++) { 
                    if(index >= index_intervel + start_index){//right
                        p.b_right.enclose(primitives[index].bbox());
                        p.right_prim_count++;
                    }
                    else{
                        p.b_left.enclose(primitives[index].bbox());
                        p.left_prim_count++;
                    }
                }
                // end of one partition
                // SAH Cost(for this axis)
                float cost = p.b_left.surface_area() / box.surface_area() * (float)p.left_prim_count +
                             p.b_right.surface_area() / box.surface_area() * (float)p.right_prim_count;
                if(cost < best_cost_ofxyz[axis]) {
                    best_cost_ofxyz[axis] = cost;
                    best_partition_thisaxis = p;
                }
            }

            best_partition[axis] = best_partition_thisaxis;
            //TO DO primitives need to sort() based on partition(
            // std::partition(primitives.begin(),primitives.end(),[](){ return bbox().center() < })
        }
        best_cost = std::min(best_cost_ofxyz[0], best_cost_ofxyz[1], best_cost_ofxyz[2]);
        int axis = (best_cost == best_cost_ofxyz[0])?(0):(   (best_cost == best_cost_ofxyz[1])? (1) : (2)  );
            nodes[root_idx].l = (size_t)nodes.size();
            nodes[root_idx].r = (size_t)nodes.size() + 1;
            
            //need to sort primitive to make sure
            
            // left
            new_node(best_partition[axis].b_left, nodes[root_idx].start, best_partition[axis].left_prim_count, 0, 0);
            // right
            new_node(best_partition[axis].b_right, nodes[root_idx].start + best_partition[axis].left_prim_count, best_partition[axis].right_prim_count,0, 0);
            

        cur_idx++;
    }

}

// template<typename Primitive> //, Vec2& times
// void BVH<Primitive>::find_closest_hit(const Ray& ray, Node& node, Trace& ret) const{
//     //if node is leaf
//     if(node.is_leaf){
//         //Each primitive in this node
//         //size is how many primitives under this node(include sub tree)
//         for(int index = node.start;index < node.size;index++) {
//             Trace hit = primitives[index].hit(ray);
//             ret = Trace::min(ret, hit);//return the mini Trace
//         }
//     }
//     else{
//         //ray.dist_bounds.x

//         //need a ret range
//         Vec2 times1 = ray.dist_bounds/(ray.dir.norm() );//std::abs((t2 * ray.dir).norm());
//         Vec2 times2 = times1;
//         bool hit_left = nodes[node.l].bbox.hit(ray, times1);
//         bool hit_right = nodes[node.r].bbox.hit(ray, times2);
//         //if left child bbox is closer thant Right, it doesn't mean that left has the closest primitive
//         if(!hit_left && !hit_right)continue;
//         else if(hit_left && hit_right){
            
//             float closer_index = (times1.x < times2.x) ? node.l : node.r;
//             float second_index = (times1.x < times2.x) ? node.r : node.l;
//             float second_hit_time = (times1.x < times2.x) ? times2.x : times1.x;
//             float second_hit_distance = std::abs((second_hit_time * ray.dir).norm());
            
//             find_closest_hit(ray, nodes[closer_index], ret); 
//             if(second_hit_distance < ray.dist_bounds.x){
//                 find_closest_hit(ray, nodes[second_index], ret); 
//             }

//         }
            
            


//         }
// }


template<typename Primitive> Trace BVH<Primitive>::hit(const Ray& ray) const {

    // TODO (PathTracer): Task 3
    // Implement ray - BVH intersection test. A ray intersects
    // with a BVH aggregate if and only if it intersects a primitive in
    // the BVH that is not an aggregate.

    // The starter code simply iterates through all the primitives.
    // Again, remember you can use hit() on any Primitive value.
    
    Trace ret;
    if(nodes.empty())return ret;
    for(const Primitive& prim : primitives) {
        Trace hit = prim.hit(ray);
        ret = Trace::min(ret, hit);//return the mini Trace
    }
    

    
    //find_closest_hit(ray, nodes[0], ret);
    
    
    
    
    
    //Level order traversal
    // std::queue<Node> q;
    // q.push(nodes[0]);
    // Vec2 times(0,FLT_MAX);
    // while(!q.empty()){
    //     auto node = q.front();
    //     q.pop();

    // }
    // for(const Node& node : nodes) {
        


        
    // }




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
