
#include "../rays/bvh.h"
#include "debug.h"
#include <stack>

namespace PT {


struct Partition {
    
    int axis = 0;
    BBox b_left, b_right;
    int left_prim_count = 0, right_prim_count = 0;
    float line = 0;
    // Partition(){
    //     axis = 0;
    //     b_left = BBox(); b_right= BBox();
    // }
};

// template<typename Primitive>
// bool x_cmp(const Primitive& prim1, const Primitive& prim2) {  
//     return (prim1.bbox().center().x < prim2.bbox().center().x);
// }
// template<typename Primitive>
// bool y_cmp(const Primitive& prim1, const Primitive& prim2) {  
//     return (prim1.bbox().center().y < prim2.bbox().center().y);
// }
// template<typename Primitive>
// bool z_cmp(const Primitive& prim1, const Primitive& prim2) {  
//     return (prim1.bbox().center().z < prim2.bbox().center().z);
// }

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
    while(true) {//++cur_idx
        if(cur_idx >= (int)nodes.size()) return;
        if(nodes[cur_idx].size <= max_leaf_size) {
            cur_idx++;
            continue;
        }
        //primitive start index and end index
        int start_index = (int)nodes[cur_idx].start;
        int end_index = (int)nodes[cur_idx].start + (int)nodes[cur_idx].size;
        // For axis x,y,z
        float best_cost = FLT_MAX;
        Vec3 best_cost_ofxyz(FLT_MAX, FLT_MAX, FLT_MAX);
        std::vector<Partition> best_partition(3); // self-define
        for(int axis = 0; axis < 3; axis++) {
            Partition best_partition_thisaxis;
            // Sort the primitives in ascending order, override the compare using lambda
            //based on this sub list
            

            // Initialize buckets

            float bucket_interval = (nodes[cur_idx].bbox.max[axis] - nodes[cur_idx].bbox.min[axis]) / (float)bucket_number;
            for(float median = nodes[cur_idx].bbox.min[axis] + bucket_interval; median < nodes[cur_idx].bbox.max[axis]; median += bucket_interval) { // the interval is 1
                // range: bbox.min[axis] -> line; line -> bbox.max[axis];
                // from axis = line to (line+1 or box.max[axis](when out of range))
                // it : the first one of the right part

                auto it = std::partition(primitives.begin() + start_index, primitives.begin() + end_index, [median, axis](auto& a){//const Primitive &  ,&axis
                    float center = a.bbox().center()[axis];
                    return center < median;
                });
                int median_index = (int)(it - (primitives.begin() ));//first index of right part + start_index 

                //index_intervel could be primitive.end() out of range
                Partition p; p.line = median;
                for(int index = start_index; index < end_index; index++) { 
                    if(index >= median_index){//right 
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
                float cost = p.b_left.surface_area() / nodes[cur_idx].bbox.surface_area() * (float)p.left_prim_count +
                             p.b_right.surface_area() / nodes[cur_idx].bbox.surface_area() * (float)p.right_prim_count
                             + 1.0f;
                if(cost < best_cost_ofxyz[axis]) {
                    best_cost_ofxyz[axis] = cost;
                    best_partition_thisaxis = p;
                }
            }
            //best_partition.insert( best_partition.begin() + axis, best_partition_thisaxis);
            best_partition[axis] = best_partition_thisaxis;
            //TO DO primitives need to sort() based on partition(
            // std::partition(primitives.begin(),primitives.end(),[](){ return bbox().center() < })
        }
        best_cost = std::min(best_cost_ofxyz[0], std::min(best_cost_ofxyz[1], best_cost_ofxyz[2]));
        int bestaxis = (best_cost == best_cost_ofxyz[0])?(0):(   (best_cost == best_cost_ofxyz[1])? (1) : (2)  );
        nodes[cur_idx].l = (size_t)nodes.size();//left child index in nodes
        nodes[cur_idx].r = (size_t)nodes.size() + 1;//right child index in nodes
        float line = best_partition[bestaxis].line;
            //need to sort primitive to make sure
            std::partition(primitives.begin() + start_index, primitives.begin() + end_index, [bestaxis, line](auto& a){// auto best_it = //const Primitive &  ,&axis
                    float center = a.bbox().center()[bestaxis];
                    return center < line;
            });
            // left
            //bbox, start index in primitives, size_t size, size_t l, size_t r
        new_node(best_partition[bestaxis].b_left, nodes[cur_idx].start, best_partition[bestaxis].left_prim_count, 0, 0);
        //new_node(best_partition[bestaxis].b_left, nodes[cur_idx].start, (int)(best_it - (primitives.begin() + start_index)), 0, 0);
            // right
        new_node(best_partition[bestaxis].b_right, nodes[cur_idx].start + best_partition[bestaxis].left_prim_count, best_partition[bestaxis].right_prim_count,0, 0);
            

        cur_idx++;
    }

}

template<typename Primitive>
Trace BVH<Primitive>::find_closest_hit(const Ray& ray, size_t root, Vec2 &times) const{
    Trace ret;
    //if node is leaf
    if(nodes[root].is_leaf()){
        //Each primitive in this node
        //size is how many primitives under this node(include sub tree)
        for(size_t index = (size_t)nodes[root].start;index < (size_t)nodes[root].start + (size_t)nodes[root].size;index++) {
            Trace hit = primitives[index].hit(ray);
            ret = Trace::min(ret, hit);//return the mini Trace
        }
        times.x = ret.distance;
        return ret;
    }
    else{
        //need a ret range
        Vec2 times1 = times;//ray.dist_bounds/(ray.dir.norm() );//std::abs((t2 * ray.dir).norm());
        Vec2 times2 = times;
        bool hit_left = nodes[nodes[root].l].bbox.hit(ray, times1);
        bool hit_right = nodes[nodes[root].r].bbox.hit(ray, times2);
        //if left child bbox is closer thant Right, it doesn't mean that left has the closest primitive
        if(!hit_left && !hit_right) return ret;
        else{
            int closer_index, second_index;
            bool hitboth = false;
            //float first_hit_time, second_hit_time;//, first_hit_distance, second_hit_distance;
            Vec2 cur_close_t = ray.dist_bounds;
            Vec2 cur_far_t = ray.dist_bounds;
            if(hit_left && hit_right){
                hitboth = true;
                if(times1.x < times2.x){
                    closer_index = (int)nodes[root].l; second_index = (int)nodes[root].r; //first_hit_time = times1.x; second_hit_time = times2.x;
                    cur_close_t = times1;cur_far_t = times2;
                }
                else{
                    closer_index = (int)nodes[root].r; second_index = (int)nodes[root].l; //first_hit_time = times2.x; second_hit_time = times1.x;
                    cur_close_t = times2;cur_far_t = times1;
                }
            }
            else if (hit_left){
                closer_index = (int)nodes[root].l; second_index = (int)nodes[root].r; //first_hit_time = times1.x; second_hit_time = times2.y;
                cur_close_t = times1;
            }
            else{
                closer_index = (int)nodes[root].r; second_index = (int)nodes[root].l; //first_hit_time = times2.x; second_hit_time = times1.y;
                cur_close_t = times2;
            }
            //second_hit_distance = std::abs((second_hit_time * ray.dir).norm());
            //first_hit_distance = std::abs(( first_hit_time * ray.dir).norm());
            //ray.dist_bounds = cur_close_t;//first_hit_distance;
            ret = find_closest_hit(ray, closer_index, cur_close_t); 
            if(cur_far_t.x < ray.dist_bounds.y && hitboth){ // cur_far_t.x < ray.dist_bounds.y
                Trace hit  = find_closest_hit(ray, second_index, cur_far_t); 
                ret = Trace::min(ret,hit);
            }
            return ret;
        }
    }
}


template<typename Primitive>
Trace BVH<Primitive>::hit(const Ray& ray) const {//const;

    // TODO (PathTracer): Task 3
    // Implement ray - BVH intersection test. A ray intersects
    // with a BVH aggregate if and only if it intersects a primitive in
    // the BVH that is not an aggregate.

    // The starter code simply iterates through all the primitives.
    // Again, remember you can use hit() on any Primitive value.
    
    Trace ret;
    if(nodes.empty())return ret;

    // for(const Primitive& prim : primitives) {
    //     Trace hit = prim.hit(ray);
    //     ret = Trace::min(ret, hit);
    // }
    
    size_t root = 0;
    Vec2 time_initial = ray.dist_bounds/(ray.dir.norm() );
    //Node *node = nodes[0];
    Trace result;
    //result = find_closest_hit(ray, root, time_initial, ret);
    result = find_closest_hit(ray, root, time_initial);
    
    
    
    
    
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




    return result;
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
