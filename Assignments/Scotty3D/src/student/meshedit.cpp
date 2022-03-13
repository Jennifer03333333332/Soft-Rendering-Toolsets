
#include <queue>
#include <set>
#include <unordered_map>

#include "../geometry/halfedge.h"
#include "debug.h"

/******************************************************************
*********************** Local Operations **************************
******************************************************************/

/* Note on local operation return types:

    The local operations all return a std::optional<T> type. This is used so that your
    implementation can signify that it does not want to perform the operation for
    whatever reason (e.g. you don't want to allow the user to erase the last vertex).

    An optional can have two values: std::nullopt, or a value of the type it is
    parameterized on. In this way, it's similar to a pointer, but has two advantages:
    the value it holds need not be allocated elsewhere, and it provides an API that
    forces the user to check if it is null before using the value.

    In your implementation, if you have successfully performed the operation, you can
    simply return the required reference:

            ... collapse the edge ...
            return collapsed_vertex_ref;

    And if you wish to deny the operation, you can return the null optional:

            return std::nullopt;

    Note that the stubs below all reject their duties by returning the null optional.
*/

/*
    This method splits the given edge in half, but does not split the
    adjacent faces. Returns an iterator to the new vertex which splits
    the original edge.

    Example function for how to go about implementing local operations
*/
std::optional<Halfedge_Mesh::VertexRef> Halfedge_Mesh::bisect_edge(EdgeRef e) {

    // Phase 1: collect all elements
    HalfedgeRef h = (e->halfedge()->is_boundary()) ? e->halfedge()->twin() : e->halfedge();
    HalfedgeRef ht = h->twin();
    HalfedgeRef preh = h;
    HalfedgeRef nexht = ht->next();
    do {
        preh = preh->next();
    } while(preh->next() != h);
    Vec3 vpos = (h->vertex()->pos + ht->vertex()->pos) / 2;

    // Phase 2: Allocate new elements
    VertexRef c = new_vertex();
    c->pos = vpos;
    HalfedgeRef hn = new_halfedge();
    HalfedgeRef hnt = new_halfedge();
    EdgeRef e0 = new_edge();

    // The following elements aren't necessary for the bisect_edge, but they are here to demonstrate
    // phase 4
    FaceRef f_not_used = new_face();
    HalfedgeRef h_not_used = new_halfedge();

    // Phase 3: Reassign elements
    e0->halfedge() = hn;
    hn->twin() = hnt;
    hn->edge() = e0;
    hn->vertex() = h->vertex();
    hn->face() = h->face();
    preh->next() = hn;
    hn->next() = h;
    h->vertex() = c;
    ht->next() = hnt;
    c->halfedge() = h;
    hn->vertex()->halfedge() = hn;
    c->is_new = true;

    // example of set_neighbors:
    // condenses hnt->next() = nexht; hnt->twin() = hn; hnt->vertex() = c; hnt->edge() = e0;
    // hnt->face() = ht->face(); into one line
    hnt->set_neighbors(nexht, hn, c, e0, ht->face());

    // Phase 4: Delete unused elements
    erase(f_not_used);
    erase(h_not_used);

    // Phase 5: Return the correct iterator
    return c;
}

/*
    This method should replace the given vertex and all its neighboring
    edges and faces with a single face, returning the new face.
 */
std::optional<Halfedge_Mesh::FaceRef> Halfedge_Mesh::erase_vertex(Halfedge_Mesh::VertexRef v) {

    // Should I delete a vertex that would cause the face disappearing? No I think.
    if(v->on_boundary()) return v->halfedge()->face();
    // First find the halfedge neighbors

    HalfedgeRef h = v->halfedge();
    std::vector<HalfedgeRef> halfedge_array;
    std::vector<VertexRef> vertices_array;
    std::vector<FaceRef> faces_array;
    // Find all neighbors for v
    do {
        // h = h->twin()->next();
        // halfedge_array.push_back(h);

        HalfedgeRef tmp_h;
        // all hf in Each face
        for(tmp_h = h->next(); tmp_h->next() != h; tmp_h = tmp_h->next()) {
            halfedge_array.push_back(tmp_h);
            vertices_array.push_back(tmp_h->vertex());
        }
        faces_array.push_back(h->face());
        h = tmp_h->twin();
    } while(h != v->halfedge());
    // reassign
    for(int i = 0; i < (int)halfedge_array.size(); i++) {
        vertices_array[i]->halfedge() = halfedge_array[i];
        halfedge_array[i]->next() = halfedge_array[(i + 1) % halfedge_array.size()];
        halfedge_array[i]->face() = faces_array[0];
    }
    faces_array[0]->halfedge() = halfedge_array[0];

    // erase face
    for(int i = 1; i < (int)faces_array.size(); i++) {
        erase(faces_array[i]);
    }
    // find all neighbors hf of v
    h = v->halfedge();
    do {
        HalfedgeRef nh = h->twin()->next();
        erase(h->edge());
        erase(h->twin());
        erase(h);
        h = nh;
    } while(h != v->halfedge());

    // Edge case

    // error: a live hf's face was erased. So for every hf in every face, must set a new face

    erase(v);
    // std::optional<std::pair<Halfedge_Mesh::ElementRef, std::string>> test = validate();

    return faces_array[0];
}

/*
    This method should erase the given edge and return an iterator to the
    merged face.
 */
std::optional<Halfedge_Mesh::FaceRef> Halfedge_Mesh::erase_edge(Halfedge_Mesh::EdgeRef e) {

    // Can I erase boundary? No
    if(e->on_boundary()) return std::nullopt; // return e->halfedge()->face();
    // First find the halfedge
    HalfedgeRef h0 = e->halfedge();
    HalfedgeRef h3 = e->halfedge()->twin();

    // find the face
    FaceRef f1 = h3->face();
    FaceRef f0 = h0->face();
    // Edge case
    bool p1 = h0->next() == h0->twin(), p2 = h3->next() == h3->twin();
    if(f1 == f0 && (p1 || p2)) {
        if(p2) swap(h0, h3);
        HalfedgeRef last_hf0, last_hf1 = h3->next();
        for(HalfedgeRef h_iter = h0->next(); h_iter != h0; h_iter = h_iter->next()) {
            last_hf0 = h_iter;
        }

        last_hf0->next() = last_hf1;
        h0->vertex()->halfedge() = last_hf1;
        erase(e);
        erase(h3->vertex());
        erase(h0);
        erase(h3);
        return f1;
    }
    // Find the 2 vertices
    VertexRef v0 = h0->vertex();
    VertexRef v1 = h3->vertex();
    // Face
    FaceRef merged_face = new_face();
    merged_face->halfedge() = h0->next();
    v0->halfedge() = h3->next();
    v1->halfedge() = h0->next(); // what if h0 don't have next
    // Next(), vertex->halfedge()
    // All hf from the previous face should point to the new face
    HalfedgeRef last_hf0, last_hf1; //
    for(HalfedgeRef h_iter = h0->next(); h_iter != h0; h_iter = h_iter->next()) {
        // if(h_iter->next() == h_iter->twin())
        last_hf0 = h_iter;
        h_iter->face() = merged_face;
    }
    for(HalfedgeRef h_iter = h3->next(); h_iter != h3; h_iter = h_iter->next()) {
        last_hf1 = h_iter;
        h_iter->face() = merged_face;
    }
    last_hf0->next() = h3->next();
    last_hf1->next() = h0->next();

    // a vertex's hf is erased

    // erase faces and hfs and this edge
    // if after erase, the vertice don't have edge any more, it should be erase
    erase(f1);
    erase(f0);
    erase(e);

    erase(h0);
    erase(h3); // what if there are some vertex contained this hf?

    if(v0->degree() < 1) erase(v0);
    if(v1->degree() < 1) erase(v1);
    // std::optional<std::pair<Halfedge_Mesh::ElementRef, std::string>> test = validate();

    return merged_face;
}

/*
    This method should collapse the given edge and return an iterator to
    the new vertex created by the collapse.
*/
std::optional<Halfedge_Mesh::VertexRef> Halfedge_Mesh::collapse_edge(Halfedge_Mesh::EdgeRef e) {

    // collapse boundaries is ok
    // Edge case : after collapse, no face left in this mesh / this mesh start with 2 faces / no
    // collapse if other edges are all on boundary
    if(e->on_boundary()) return std::nullopt;
    // without h0, h3
    std::vector<HalfedgeRef> hf_inFace0, hf_inFace1; //, hf_outside;
    // std::vector<VertexRef> vertex_array;
    // std::vector<EdgeRef> edges_array;

    // First find the halfedge
    HalfedgeRef h0 = e->halfedge();
    HalfedgeRef h3 = e->halfedge()->twin();

    // find the face
    FaceRef f1 = h3->face();
    FaceRef f0 = h0->face();

    if(f0->degree() <= 3 &&
       f0->is_boundary()) { // tetrahedron's face won't be boundary(manifold assumption)
        return std::nullopt;
    }
    if(f1->degree() <= 3 && f1->is_boundary()) {
        return std::nullopt;
    }
    // int test = (int)faces.size();
    // Find the 2 vertices
    VertexRef v0 = h0->vertex();
    VertexRef v1 = h3->vertex();
    VertexRef v = new_vertex();
    v->pos = (v0->pos + v1->pos) / 2;
    // deal with next
    // Store 2 arrays for face0 and face1
    for(HalfedgeRef h_iter = h3->next(); h_iter != h3; h_iter = h_iter->next()) {
        // h3_last = h_iter;
        hf_inFace1.push_back(h_iter);
    }
    for(HalfedgeRef h_iter = h0->next(); h_iter != h0; h_iter = h_iter->next()) {
        // h0_last = h_iter;
        hf_inFace0.push_back(h_iter);
    }
    // get all the edges/hfs with v0 and v1, change it all from v0/v1 to v
    {
        // change//std::vector<HalfedgeRef> hf_withv0,hf_withv1;
        HalfedgeRef h = v0->halfedge();
        do {
            // h = h->twin();
            // hf_withv0.push_back(h);
            h->vertex() = v;
            h = h->twin()->next();
            // hf_withv0.push_back(h);//start point = v0
        } while(h != v0->halfedge());
        h = v1->halfedge();
        do {
            h->vertex() = v;
            h = h->twin()->next();
            // hf_withv1.push_back(h);//start point = v1
        } while(h != v1->halfedge());
    }
    // 1 erase edge e //erase_edge(e);
    // deal with next
    hf_inFace1[hf_inFace1.size() - 1]->next() = h3->next();
    hf_inFace0[hf_inFace0.size() - 1]->next() = h0->next();
    // before erase h0,h3
    v0->halfedge() = hf_inFace0[hf_inFace0.size() - 1]->twin();
    v1->halfedge() = hf_inFace1[hf_inFace1.size() - 1]->twin();
    erase(h0);
    erase(h3);
    erase(e);

    // if there are 2-degree face, return null or replace it with line
    erase(v0);
    erase(v1);

    // stay face or collapse face
    //   f0
    if(hf_inFace0.size() + 1 >
       3) { // hf_inFace0.size() + 1 = edges number of f0// if it's not triangle
        f0->halfedge() = hf_inFace0[0];
        // bind next() in face
        for(int i = 0; i < (int)hf_inFace0.size(); i++) {
            if(i == 0)
                hf_inFace0[hf_inFace0.size() - 1]->next() = hf_inFace0[0];
            else {
                hf_inFace0[i - 1]->next() = hf_inFace0[i];
            }
        }
        v->halfedge() = hf_inFace0[0];
    } else { // triangle
        HalfedgeRef a = hf_inFace0[0], b = hf_inFace0[1], a_twin = a->twin(), b_twin = b->twin();
        VertexRef tmp_v0 = a_twin->vertex();
        VertexRef tmp_v1 = b_twin->vertex();
        // bind twin(). keep atwin and btwin, erase a ,b
        a_twin->twin() = b_twin;
        b_twin->twin() = a_twin;
        // before erase(b->edge());
        b_twin->edge() = a_twin->edge();
        // before erase a: twin,next,vertex,edge,face
        a->vertex()->halfedge() = b_twin;
        a->edge()->halfedge() = a_twin;
        a_twin->vertex()->halfedge() = a_twin;
        // before erase b :
        b->vertex()->halfedge() = b_twin;
        b_twin->vertex()->halfedge() = b_twin;

        tmp_v0->halfedge() = a_twin;
        tmp_v1->halfedge() = b_twin;

        erase(b->face());
        erase(b->edge());
        erase(a);
        erase(b);
        erase(f0);
    }

    // f1
    if(hf_inFace1.size() + 1 >
       3) { // hf_inFace0.size() + 1 = edges number of f0// if it's not triangle
        f1->halfedge() = hf_inFace1[0];
        // bind next() in face
        for(int i = 0; i < (int)hf_inFace1.size(); i++) {
            if(i == 0)
                hf_inFace1[hf_inFace1.size() - 1]->next() = hf_inFace1[0];
            else {
                hf_inFace1[i - 1]->next() = hf_inFace1[i];
            }
        }
    } else { // triangle
        HalfedgeRef a = hf_inFace1[0], b = hf_inFace1[1], a_twin = a->twin(), b_twin = b->twin();
        VertexRef tmp_v0 = a_twin->vertex();
        VertexRef tmp_v1 = b_twin->vertex();
        // bind twin(). keep atwin and btwin, erase a ,b
        a_twin->twin() = b_twin;
        b_twin->twin() = a_twin;
        // before erase(b->edge());
        b_twin->edge() = a_twin->edge();
        // before erase a: twin,next,vertex,edge,face
        a->vertex()->halfedge() = b_twin;
        a->edge()->halfedge() = a_twin;
        a_twin->vertex()->halfedge() = a_twin;
        // before erase b :
        b->vertex()->halfedge() = b_twin;
        b_twin->vertex()->halfedge() = b_twin;

        tmp_v0->halfedge() = a_twin;
        tmp_v1->halfedge() = b_twin;

        erase(b->face());
        erase(b->edge());
        erase(a);
        erase(b);
        erase(f1);
    }

    // I erase() f? so maybe I can use the face number? will

    // what to check if there are edge duplicate
    //  std::set<std::pair<unsigned int, unsigned int>> edge_ids;
    //  for (EdgeRef e_iter = edges.begin(); e_iter != edges.end(); e_iter++) {
    //      unsigned int l = e_iter->halfedge()->vertex()->id();
    //      unsigned int r = e_iter->halfedge()->twin()->vertex()->id();
    //      if(l == r) {
    //          return std::nullopt;
    //      }
    //      auto entry = edge_ids.find({l, r});
    //      if(entry != edge_ids.end()) {
    //          return std::nullopt;
    //      }
    //      edge_ids.insert({l, r});
    //      edge_ids.insert({r, l});
    //  }
    // std::optional<std::pair<Halfedge_Mesh::ElementRef, std::string>> test = validate();

    return v;
}

/*
    This method should collapse the given face and return an iterator to
    the new vertex created by the collapse.
*/
std::optional<Halfedge_Mesh::VertexRef> Halfedge_Mesh::collapse_face(Halfedge_Mesh::FaceRef f) {

    (void)f;
    return std::nullopt;
}

/*
    This method should flip the given edge counter-clockwise and return an iterator to the
    flipped edge.
*/
std::optional<Halfedge_Mesh::EdgeRef> Halfedge_Mesh::flip_edge(Halfedge_Mesh::EdgeRef e) {

    // Task 1
    // Edge case : tetrahedron, e not on boundary??

    if(e->on_boundary()) { // can't flip the edge on the boundary
        return std::nullopt;
    }
    std::vector<HalfedgeRef> hf_inFace0, hf_inFace1, hf_outside;
    // First find the halfedge
    HalfedgeRef h0 = e->halfedge();
    HalfedgeRef h3 = e->halfedge()->twin();

    // find the face
    FaceRef f1 = h3->face();
    FaceRef f0 = h0->face();
    // Find the 2 vertices
    //  VertexRef v0 =  h0->vertex();
    //  VertexRef v1 =  h3->vertex();
    // iterate hf in face 1, 0
    HalfedgeRef hf_iter = h3;
    do {
        hf_inFace1.push_back(hf_iter);
        hf_iter = hf_iter->next();
    } while(hf_iter != h3);

    hf_iter = h0;
    do {
        hf_inFace0.push_back(hf_iter);
        hf_iter = hf_iter->next();
    } while(hf_iter != h0);

    // what's in vertex_array : v0,v2,v1,v3(triangle mesh)
    // how to guide it's ccw: go next() in one face

    // Our new edge's vertex:v2,v3(h3->next()->next()->vertex(),h0->next()->next()->vertex())
    // what changed: v2,v3's hf, e's vertex, h0,h3's vertex, f0,f1's vertex
    // reassign
    // set_neighbors(HalfedgeRef next, HalfedgeRef twin, VertexRef vertex, EdgeRef edge,FaceRef
    // face)
    h3->set_neighbors(hf_inFace1[2], h0, hf_inFace0[2]->vertex(), e, f1);
    h0->set_neighbors(hf_inFace0[2], h3, hf_inFace1[2]->vertex(), e, f0);
    f1->halfedge() = h3;
    f0->halfedge() = h0;

    // Adjust face array
    HalfedgeRef tmp = hf_inFace0[1];
    hf_inFace0.erase(hf_inFace0.begin() + 1);
    hf_inFace0.push_back(hf_inFace1[1]);

    hf_inFace1.erase(hf_inFace1.begin() + 1);
    hf_inFace1.push_back(tmp);
    // Iterate based on hf array, fix next()
    for(int iter = 0; iter < (int)hf_inFace0.size(); iter++) {
        if(iter > 0)
            hf_inFace0[iter - 1]->next() = hf_inFace0[iter];
        else {
            hf_inFace0[hf_inFace0.size() - 1]->next() = hf_inFace0[0];
        }
        hf_inFace0[iter]->face() = f0;
    }
    for(int iter = 0; iter < (int)hf_inFace1.size(); iter++) {
        if(iter > 0)
            hf_inFace1[iter - 1]->next() = hf_inFace1[iter];
        else {
            hf_inFace1[hf_inFace1.size() - 1]->next() = hf_inFace1[0];
        }
        hf_inFace1[iter]->face() = f1;
    }
    // what if after flip it duplicate with other edge?
    // hf_inFace0[2]->next() = hf_inFace1[1];
    // error how do you know this face don't have next hf?
    // v0,v1 is not in the new face
    if(f0->degree() <= 2 && f1->degree() <= 2) return std::nullopt;
    // how to judge 3 points is on the same edge: using the same vertex
    // how to judge if there are edge duplicate
    std::set<std::pair<unsigned int, unsigned int>> edge_ids;
    for(EdgeRef e_iter = edges.begin(); e_iter != edges.end(); e_iter++) {
        unsigned int l = e_iter->halfedge()->vertex()->id();
        unsigned int r = e_iter->halfedge()->twin()->vertex()->id();
        if(l == r) {
            return std::nullopt;
        }
        auto entry = edge_ids.find({l, r});
        if(entry != edge_ids.end()) {
            return std::nullopt;
        }
        edge_ids.insert({l, r});
        edge_ids.insert({r, l});
    }
    // std::optional<std::pair<Halfedge_Mesh::ElementRef, std::string>> test = validate();
    return e;
}

/*
    This method should split the given edge and return an iterator to the
    newly inserted vertex. The halfedge of this vertex should point along
    the edge that was split, rather than the new edges.
*/
std::optional<Halfedge_Mesh::VertexRef> Halfedge_Mesh::split_edge(Halfedge_Mesh::EdgeRef e) {

    if(e->on_boundary()) { // can't flip the edge on the boundary
        return std::nullopt;
    }

    std::vector<HalfedgeRef> hf_inFace0, hf_inFace1;
    std::vector<HalfedgeRef> newhf_f0, newhf_f1, newhf_f2, newhf_f3;
    // First find the halfedge
    HalfedgeRef h0 = e->halfedge();
    HalfedgeRef h3 = e->halfedge()->twin();
    EdgeRef olde0 = h0->edge();
    EdgeRef olde3 = h3->edge();
    // find the face
    FaceRef f1 = h3->face();
    FaceRef f0 = h0->face();
    // Find the 2 vertices
    VertexRef v0 = h0->vertex();
    VertexRef v1 = h3->vertex();

    // Start the reassign

    // For the 2 new edges and For the original edge that was split
    EdgeRef new_e0 = new_edge(), new_e1 = new_edge();
    HalfedgeRef new_hf1 = new_halfedge(), new_hf1_twin = new_halfedge(), new_hf0 = new_halfedge(),
                new_hf0_twin = new_halfedge();

    EdgeRef new_e2 = new_edge(), new_e3 = new_edge();
    HalfedgeRef new_hf2 = new_halfedge(), new_hf2_twin = new_halfedge(), new_hf3 = new_halfedge(),
                new_hf3_twin = new_halfedge();
    FaceRef f2 = new_face(), f3 = new_face();

    // iterate hf in face 1
    hf_inFace1.push_back(h3);
    for(HalfedgeRef hf_iter = h3->next(); hf_iter != h3; hf_iter = hf_iter->next()) {
        hf_inFace1.push_back(hf_iter);
    }
    // iterate hf in face 0
    hf_inFace0.push_back(h0);
    for(HalfedgeRef hf_iter = h0->next(); hf_iter != h0; hf_iter = hf_iter->next()) {
        hf_inFace0.push_back(hf_iter);
    }
    // for triangle only
    if(hf_inFace1.size() != 3 && hf_inFace0.size() != 3) return std::nullopt;
    VertexRef v3 = hf_inFace0[2]->vertex();
    VertexRef v2 = hf_inFace1[2]->vertex();
    // Storage ends
    //  new vertex
    VertexRef new_v = new_vertex();
    new_v->pos = (v3->pos + v2->pos) / 2;

    // assign hf to edge
    new_e2->halfedge() = new_hf2;
    new_e1->halfedge() = new_hf1;
    new_e0->halfedge() = new_hf0;
    new_e3->halfedge() = new_hf3;
    // assign hf to face
    f0->halfedge() = new_hf0;
    f1->halfedge() = new_hf1;
    f2->halfedge() = new_hf2;
    f3->halfedge() = new_hf3;
    // assign hf to vertex
    v0->halfedge() = new_hf3;
    v1->halfedge() = new_hf2;
    v2->halfedge() = new_hf1;
    v3->halfedge() = new_hf0;
    new_v->halfedge() = new_hf2_twin; //! forget the hf of new vertex
    // assign everything to hf
    new_hf0->set_neighbors(new_hf2_twin, new_hf0_twin, v3, new_e0, f0);
    new_hf0_twin->set_neighbors(hf_inFace0[2], new_hf0, new_v, new_e0, f3);
    new_hf1->set_neighbors(new_hf3_twin, new_hf1_twin, v2, new_e1, f1);
    new_hf1_twin->set_neighbors(hf_inFace1[2], new_hf1, new_v, new_e1, f2);
    new_hf2->set_neighbors(new_hf1_twin, new_hf2_twin, v1, new_e2, f2);
    new_hf2_twin->set_neighbors(hf_inFace0[1], new_hf2, new_v, new_e2, f0);
    new_hf3->set_neighbors(new_hf0_twin, new_hf3_twin, v0, new_e3, f3);
    new_hf3_twin->set_neighbors(hf_inFace1[1], new_hf3, new_v, new_e3, f1);
    //! forget the face of olf hf
    hf_inFace0[1]->face() = f0;
    hf_inFace0[2]->face() = f3;
    hf_inFace1[1]->face() = f1;
    hf_inFace1[2]->face() = f2;
    // next() for the old hf
    hf_inFace0[1]->next() = new_hf0;
    hf_inFace0[2]->next() = new_hf3;
    hf_inFace1[1]->next() = new_hf1;
    hf_inFace1[2]->next() = new_hf2;

    // erase the old edges
    erase(olde0);
    erase(olde3);
    // erase hf: h0,h3
    erase(h0);
    erase(h3);
    // erase face: in this case no need
    // vertex: nothing changed
    // std::optional<std::pair<Halfedge_Mesh::ElementRef, std::string>> test = validate();
    return new_v;
}

/*
    This method should insets a vertex into the given face, returning a pointer to the new center
   vertex
*/
std::optional<Halfedge_Mesh::VertexRef> Halfedge_Mesh::inset_vertex(FaceRef f) {
    (void)f;
    return std::nullopt;
}

/*
    This method should inset a face into the given face, returning a pointer to the new face.
*/
std::optional<Halfedge_Mesh::FaceRef> Halfedge_Mesh::inset_face(Halfedge_Mesh::FaceRef f) {

    // hint: use bevel_face positions as a helper function here
    (void)f;
    return std::nullopt;
}

/*
    This method should bevel a vertex and inserts a vertex into the new vertex, returning a pointer
   to that vertex
*/
std::optional<Halfedge_Mesh::VertexRef> Halfedge_Mesh::extrude_vertex(VertexRef v) {
    (void)v;
    return std::nullopt;
}

/* Note on the beveling process:

    Each of the bevel_vertex, bevel_edge, and bevel_face functions do not represent
    a full bevel operation. Instead, they should update the _connectivity_ of
    the mesh, _not_ the positions of newly created vertices. In fact, you should set
    the positions of new vertices to be exactly the same as wherever they "started from."

    When you click on a mesh element while in bevel mode, one of those three functions
    is called. But, because you may then adjust the distance/offset of the newly
    beveled face, we need another method of updating the positions of the new vertices.

    This is where bevel_vertex_positions, bevel_edge_positions, and
    bevel_face_positions come in: these functions are called repeatedly as you
    move your mouse, the position of which determins the normal and tangent offset
    parameters. These functions are also passed an array of the original vertex
    positions: for bevel_vertex, it has one element, the original vertex position,
    for bevel_edge, two for the two vertices, and for bevel_face, it has the original
    position of each vertex in order starting from face->halfedge. You should use these
    positions, as well as the normal and tangent offset fields to assign positions to
    the new vertices.

    Finally, note that the normal and tangent offsets are not relative values - you
    should compute a particular new position from them, not a delta to apply.
*/

/*
    This method should replace the vertex v with a face, corresponding to
    a bevel operation. It should return the new face.  NOTE: This method is
    only responsible for updating the *connectivity* of the mesh---it does not
    need to update the vertex positions. These positions will be updated in
    Halfedge_Mesh::bevel_vertex_positions (which you also have to
    implement!)
*/
std::optional<Halfedge_Mesh::FaceRef> Halfedge_Mesh::bevel_vertex(Halfedge_Mesh::VertexRef v) {

    // Reminder: You should set the positions of new vertices (v->pos) to be exactly
    // the same as wherever they "started from."

    if(v->degree() <= 2) return std::nullopt;

    std::vector<HalfedgeRef> halfedges_array, inside;
    std::vector<VertexRef> verts;
    HalfedgeRef th = v->halfedge();
    do {
        th = th->twin();
        halfedges_array.push_back(th);
        th = th->next();
        halfedges_array.push_back(th);

        VertexRef nv = new_vertex();
        nv->pos = v->pos;
        verts.push_back(nv);
        inside.push_back(new_halfedge());
    } while(th != v->halfedge());

    FaceRef f = new_face();
    f->halfedge() = inside[0];

    int deg = v->degree();
    for(int i = 0; i < deg; i++) {
        HalfedgeRef nh = new_halfedge();
        EdgeRef ne = new_edge();

        nh->twin() = inside[i];
        inside[i]->twin() = nh;
        nh->edge() = inside[i]->edge() = ne;
        ne->halfedge() = nh;

        nh->vertex() = verts[i];
        verts[i]->halfedge() = nh;

        inside[i]->vertex() = verts[(i + 1) % deg];
        halfedges_array[i << 1 | 1]->vertex() = verts[(i + 1) % deg];

        halfedges_array[i << 1]->next() = nh;
        nh->next() = halfedges_array[i << 1 | 1];
        inside[i]->next() = inside[(i + deg - 1) % deg];

        inside[i]->face() = f;
        nh->face() = halfedges_array[i << 1]->face();
    }

    erase(v);

    return f;
}

/*
    This method should replace the edge e with a face, corresponding to a
    bevel operation. It should return the new face. NOTE: This method is
    responsible for updating the *connectivity* of the mesh only---it does not
    need to update the vertex positions. These positions will be updated in
    Halfedge_Mesh::bevel_edge_positions (which you also have to
    implement!)
*/
std::optional<Halfedge_Mesh::FaceRef> Halfedge_Mesh::bevel_edge(Halfedge_Mesh::EdgeRef e) {

    // Reminder: You should set the positions of new vertices (v->pos) to be exactly
    // the same as wherever they "started from."

    (void)e;
    return std::nullopt;
}

/*
    This method should replace the face f with an additional, inset face
    (and ring of faces around it), corresponding to a bevel operation. It
    should return the new face.  NOTE: This method is responsible for updating
    the *connectivity* of the mesh only---it does not need to update the vertex
    positions. These positions will be updated in
    Halfedge_Mesh::bevel_face_positions (which you also have to
    implement!)
*/

inline int Findh0TwinInSideFace(const int& index, const int& degree) {
    return ((index + 1) * 3 + 2) % (3 * degree); // number of new hf for each side face
}
std::optional<Halfedge_Mesh::FaceRef> Halfedge_Mesh::bevel_face(Halfedge_Mesh::FaceRef f) {

    // Reminder: You should set the positions of new vertices (v->pos) to be exactly
    // the same as wherever they "started from."

    // if (f->is_boundary()) {return f;}

    std::vector<HalfedgeRef> old_hfs, new_hfs, halfedges_array;
    std::vector<FaceRef> faces_array;      // side faces
    std::vector<VertexRef> vertices_array; // vertices in up face
    // Store hfs in old face in old_hfs
    HalfedgeRef h = f->halfedge();
    do {
        old_hfs.emplace_back(h);
        h = h->next();
    } while(h != f->halfedge());
    int degree = (int)old_hfs.size();
    // For example, 4 degree polygon -> 4 faces as side
    for(int i = 0; i < degree; i++) {

        FaceRef f0 = new_face();
        // Bind face and hf
        f0->halfedge() = old_hfs[i];
        old_hfs[i]->face() = f0;

        // bind face
        HalfedgeRef h0 = new_halfedge(), h1 = new_halfedge(), h2 = new_halfedge();

        // Store array
        // For side face: ⬆h0 →h1 ⬇h2 ⬅old_hfs[i]
        // next
        h0->next() = h1;
        h1->next() = h2;
        h2->next() = old_hfs[i];
        old_hfs[i]->next() = h0;
        // twin
        h0->face() = f0;
        h1->face() = f0;
        h2->face() = f0;
        h0->vertex() =
            old_hfs[(i + 1) % degree]->vertex(); // old_hfs[i]->next()->vertex(); old_hfs[i]->next()
                                                 // could don't have vertex

        halfedges_array.emplace_back(h0);
        halfedges_array.emplace_back(h1);
        halfedges_array.emplace_back(h2); // hf for one side face
        faces_array.emplace_back(f0);
        // For Up face:
        new_hfs.emplace_back(new_halfedge());
        vertices_array.emplace_back(new_vertex());
    }

    // assign up face
    for(int j = 0; j < degree; j++) {
        // In up face
        {
            // Bind next in up face
            new_hfs[j]->next() = new_hfs[(j + 1) % degree];
            // Bind twin with side face
            new_hfs[j]->twin() = halfedges_array[j * 3 + 1]; // h1
            halfedges_array[j * 3 + 1]->twin() = new_hfs[j];
            // Bind vertex in up face
            new_hfs[j]->vertex() = vertices_array[j];
            vertices_array[j]->halfedge() = new_hfs[j];
            vertices_array[j]->pos = old_hfs[j]->vertex()->pos; // set as the start position
            // Bind edge in up face
            new_hfs[j]->edge() = halfedges_array[j * 3 + 1]->edge() = new_edge();
            new_hfs[j]->edge()->halfedge() = new_hfs[j];
            // Bind face
            new_hfs[j]->face() = f;
            f->halfedge() = new_hfs[j];
        }
        // In side faces
        {
            // next(): done
            // set twin(!!!)
            halfedges_array[j * 3]->twin() = halfedges_array[Findh0TwinInSideFace(j, degree)];
            halfedges_array[Findh0TwinInSideFace(j, degree)]->twin() = halfedges_array[j * 3];

            // Bind vertex in side face
            // Vertex
            halfedges_array[j * 3 + 1]->vertex() = vertices_array[(j + 1) % old_hfs.size()];
            halfedges_array[j * 3 + 2]->vertex() = vertices_array[j];
            // For side edges ( vertical )

            // Bind edge
            halfedges_array[j * 3]->edge() =
                halfedges_array[Findh0TwinInSideFace(j, degree)]->edge() = new_edge();
            halfedges_array[j * 3]->edge()->halfedge() = halfedges_array[j * 3];
            // Bind face: Done
        }
    }

    // std::optional<std::pair<Halfedge_Mesh::ElementRef, std::string>> test = validate();
    return f;
}

/*
    Compute new vertex positions for the vertices of the beveled vertex.

    These vertices can be accessed via new_halfedges[i]->vertex()->pos for
    i = 1, ..., new_halfedges.size()-1.

    The basic strategy here is to loop over the list of outgoing halfedges,
    and use the original vertex position and its associated outgoing edge
    to compute a new vertex position along the outgoing edge.
*/
void Halfedge_Mesh::bevel_vertex_positions(const std::vector<Vec3>& start_positions,
                                           Halfedge_Mesh::FaceRef face, float tangent_offset) {

    std::vector<HalfedgeRef> new_halfedges;
    auto h = face->halfedge();
    do {
        new_halfedges.emplace_back(h);
        h = h->next();
    } while(h != face->halfedge());

    // (void)new_halfedges;
    // (void)start_positions;
    // (void)face;
    // (void)tangent_offset;
    tangent_offset = (-1.0f) * tangent_offset; // left and right was reverse
    // Store new_halfedges
    std::vector<Vec3> before_pos;
    for(auto e : new_halfedges) {
        before_pos.emplace_back(
            e->twin()->next()->next()->vertex()->pos); // vertex on the old edge, guide the
                                                       // direction
    }

    int degree = (int)face->degree();
    std::vector<Vec3> delta; // delta of vector's / change
    for(int i = 0; i < degree; i++) {
        Vec3 dv = before_pos[i] - start_positions[i];
        if(dv.norm_squared() != 0) dv.normalize();
        dv *= tangent_offset;
        // if get larger, Distance between(current and original edge) can not > original edge's
        // length
        if(tangent_offset > 0 &&
           (new_halfedges[i]->vertex()->pos - before_pos[i]).norm_squared() <= dv.norm_squared())
            return;
        // if get smaller, Distance between(current and original edge) can not cause reverse
        if(tangent_offset < 0 &&
           (new_halfedges[i]->vertex()->pos - start_positions[i]).norm_squared() <=
               dv.norm_squared())
            return;
        delta.emplace_back(dv);
    }

    for(int i = 0; i < degree; i++) {
        new_halfedges[i]->vertex()->pos += delta[i];
    }
}

/*
    Compute new vertex positions for the vertices of the beveled edge.

    These vertices can be accessed via new_halfedges[i]->vertex()->pos for
    i = 1, ..., new_halfedges.size()-1.

    The basic strategy here is to loop over the list of outgoing halfedges,
    and use the preceding and next vertex position from the original mesh
    (in the orig array) to compute an offset vertex position.

    Note that there is a 1-to-1 correspondence between halfedges in
    newHalfedges and vertex positions in start_positions. So, you can write
    loops of the form:

    for(size_t i = 0; i < new_halfedges.size(); i++)
    {
            Vector3D pi = start_positions[i]; // get the original vertex
            position corresponding to vertex i
    }
*/
void Halfedge_Mesh::bevel_edge_positions(const std::vector<Vec3>& start_positions,
                                         Halfedge_Mesh::FaceRef face, float tangent_offset) {

    std::vector<HalfedgeRef> new_halfedges;
    auto h = face->halfedge();
    do {
        new_halfedges.push_back(h);
        h = h->next();
    } while(h != face->halfedge());

    (void)new_halfedges;
    (void)start_positions;
    (void)face;
    (void)tangent_offset;
}

/*
    Compute new vertex positions for the vertices of the beveled face.

    These vertices can be accessed via new_halfedges[i]->vertex()->pos for
    i = 1, ..., new_halfedges.size()-1.

    The basic strategy here is to loop over the list of outgoing halfedges,
    and use the preceding and next vertex position from the original mesh
    (in the start_positions array) to compute an offset vertex
    position.

    Note that there is a 1-to-1 correspondence between halfedges in
    new_halfedges and vertex positions in start_positions. So, you can write
    loops of the form:

    for(size_t i = 0; i < new_halfedges.size(); i++)
    {
            Vec3 pi = start_positions[i]; // get the original vertex
            position corresponding to vertex i
    }
*/
void Halfedge_Mesh::bevel_face_positions(const std::vector<Vec3>& start_positions,
                                         Halfedge_Mesh::FaceRef face, float tangent_offset,
                                         float normal_offset) {

    if(tangent_offset == 0 && normal_offset == 0) return;
    tangent_offset = (-1.0f) * tangent_offset; // left and right was reverse
    if(flip_orientation) normal_offset = -normal_offset;
    std::vector<HalfedgeRef> new_halfedges;
    auto h = face->halfedge();
    do {
        new_halfedges.push_back(h);
        h = h->next();
    } while(h != face->halfedge());
    // (void)new_halfedges;
    // (void)start_positions;
    // (void)face;
    // (void)tangent_offset; : delta.pos.x
    // (void)normal_offset;: delta.pos.y
    int degree = (int)face->degree();
    Vec3 center_t = face->center();
    std::vector<Vec3> before_pos; // pos of vertices
    for(auto e : new_halfedges) {
        before_pos.push_back(e->vertex()->pos);
    }
    Vec3 normal = (-1.0f) * face->normal().normalize();
    normal *= normal_offset;

    std::vector<Vec3> delta;
    for(int i = 0; i < degree; i++) {
        Vec3 dv = Vec3();
        Vec3 horizontal = (center_t - before_pos[i]).normalize() * tangent_offset;
        dv = horizontal + normal;
        // constrained for both larger and smaller
        if((new_halfedges[i]->vertex()->pos - center_t).norm_squared() <= dv.norm_squared()) return;
        delta.push_back(dv);
    }

    for(int i = 0; i < degree; i++) {
        new_halfedges[i]->vertex()->pos += delta[i];
    }
}

/*
    Updates the position of v using the given start_position
*/
void Halfedge_Mesh::extrude_vertex_position(const Vec3& start_positions,
                                            Halfedge_Mesh::FaceRef face) {
    (void)start_positions;
    (void)face;
}

/******************************************************************
*********************** Global Operations *************************
******************************************************************/

/*
    Splits all non-triangular faces into triangles.
*/
void Halfedge_Mesh::triangulate() {

    // For each face...
}

/* Note on the quad subdivision process:

        Unlike the local mesh operations (like bevel or edge flip), we will perform
        subdivision by splitting *all* faces into quads "simultaneously."  Rather
        than operating directly on the halfedge data structure (which as you've
        seen is quite difficult to maintain!) we are going to do something a bit nicer:
           1. Create a raw list of vertex positions and faces (rather than a full-
              blown halfedge mesh).
           2. Build a new halfedge mesh from these lists, replacing the old one.
        Sometimes rebuilding a data structure from scratch is simpler (and even
        more efficient) than incrementally modifying the existing one.  These steps are
        detailed below.

  Step I: Compute the vertex positions for the subdivided mesh.
        Here we're going to do something a little bit strange: since we will
        have one vertex in the subdivided mesh for each vertex, edge, and face in
        the original mesh, we can nicely store the new vertex *positions* as
        attributes on vertices, edges, and faces of the original mesh. These positions
        can then be conveniently copied into the new, subdivided mesh.
        This is what you will implement in linear_subdivide_positions() and
        catmullclark_subdivide_positions().

  Steps II-IV are provided (see Halfedge_Mesh::subdivide()), but are still detailed
  here:

  Step II: Assign a unique index (starting at 0) to each vertex, edge, and
        face in the original mesh. These indices will be the indices of the
        vertices in the new (subdivided) mesh. They do not have to be assigned
        in any particular order, so long as no index is shared by more than one
        mesh element, and the total number of indices is equal to V+E+F, i.e.,
        the total number of vertices plus edges plus faces in the original mesh.
        Basically we just need a one-to-one mapping between original mesh elements
        and subdivided mesh vertices.

  Step III: Build a list of quads in the new (subdivided) mesh, as tuples of
        the element indices defined above. In other words, each new quad should be
        of the form (i,j,k,l), where i,j,k and l are four of the indices stored on
        our original mesh elements.  Note that it is essential to get the orientation
        right here: (i,j,k,l) is not the same as (l,k,j,i).  Indices of new faces
        should circulate in the same direction as old faces (think about the right-hand
        rule).

  Step IV: Pass the list of vertices and quads to a routine that clears
        the internal data for this halfedge mesh, and builds new halfedge data from
        scratch, using the two lists.
*/

/*
    Compute new vertex positions for a mesh that splits each polygon
    into quads (by inserting a vertex at the face midpoint and each
    of the edge midpoints).  The new vertex positions will be stored
    in the members Vertex::new_pos, Edge::new_pos, and
    Face::new_pos.  The values of the positions are based on
    simple linear interpolation, e.g., the edge midpoints and face
    centroids.
*/
void Halfedge_Mesh::linear_subdivide_positions() {

    // For each vertex, assign Vertex::new_pos to
    // its original position, Vertex::pos.

    // For each edge, assign the midpoint of the two original
    // positions to Edge::new_pos.

    // For each face, assign the centroid (i.e., arithmetic mean)
    // of the original vertex positions to Face::new_pos. Note
    // that in general, NOT all faces will be triangles!
}

/*
    Compute new vertex positions for a mesh that splits each polygon
    into quads (by inserting a vertex at the face midpoint and each
    of the edge midpoints).  The new vertex positions will be stored
    in the members Vertex::new_pos, Edge::new_pos, and
    Face::new_pos. The values of the positions are based on
    the Catmull-Clark rules for subdivision.

    Note: this will only be called on meshes without boundary
*/
void Halfedge_Mesh::catmullclark_subdivide_positions() {

    // The implementation for this routine should be
    // a lot like Halfedge_Mesh:linear_subdivide_positions:(),
    // except that the calculation of the positions themsevles is
    // slightly more involved, using the Catmull-Clark subdivision
    // rules. (These rules are outlined in the Developer Manual.)

    // Faces

    // Edges

    // Vertices
}

/*
    This routine should increase the number of triangles in the mesh
    using Loop subdivision. Note: this is will only be called on triangle meshes.
*/
void Halfedge_Mesh::loop_subdivide() {

    // Each vertex and edge of the original mesh can be associated with a
    // vertex in the new (subdivided) mesh.
    // Therefore, our strategy for computing the subdivided vertex locations is to
    // *first* compute the new positions
    // using the connectivity of the original (coarse) mesh. Navigating this mesh
    // will be much easier than navigating
    // the new subdivided (fine) mesh, which has more elements to traverse.  We
    // will then assign vertex positions in
    // the new mesh based on the values we computed for the original mesh.

    // Compute new positions for all the vertices in the input mesh using
    // the Loop subdivision rule and store them in Vertex::new_pos.
    //    At this point, we also want to mark each vertex as being a vertex of the
    //    original mesh. Use Vertex::is_new for this.

    // Next, compute the subdivided vertex positions associated with edges, and
    // store them in Edge::new_pos.

    // Next, we're going to split every edge in the mesh, in any order.
    // We're also going to distinguish subdivided edges that came from splitting
    // an edge in the original mesh from new edges by setting the boolean Edge::is_new.
    // Note that in this loop, we only want to iterate over edges of the original mesh.
    // Otherwise, we'll end up splitting edges that we just split (and the
    // loop will never end!)

    // Now flip any new edge that connects an old and new vertex.

    // Finally, copy new vertex positions into the Vertex::pos.
}

/*
    Isotropic remeshing. Note that this function returns success in a similar
    manner to the local operations, except with only a boolean value.
    (e.g. you may want to return false if this is not a triangle mesh)
*/
bool Halfedge_Mesh::isotropic_remesh() {

    // Compute the mean edge length.
    // Repeat the four main steps for 5 or 6 iterations
    // -> Split edges much longer than the target length (being careful about
    //    how the loop is written!)
    // -> Collapse edges much shorter than the target length.  Here we need to
    //    be EXTRA careful about advancing the loop, because many edges may have
    //    been destroyed by a collapse (which ones?)
    // -> Now flip each edge if it improves vertex degree
    // -> Finally, apply some tangential smoothing to the vertex positions

    // Note: if you erase elements in a local operation, they will not be actually deleted
    // until do_erase or validate is called. This is to facilitate checking
    // for dangling references to elements that will be erased.
    // The rest of the codebase will automatically call validate() after each op,
    // but here simply calling collapse_edge() will not erase the elements.
    // You should use collapse_edge_erase() instead for the desired behavior.

    return false;
}

/* Helper type for quadric simplification */
struct Edge_Record {
    Edge_Record() {
    }
    Edge_Record(std::unordered_map<Halfedge_Mesh::VertexRef, Mat4>& vertex_quadrics,
                Halfedge_Mesh::EdgeRef e)
        : edge(e) {

        // Compute the combined quadric from the edge endpoints.
        // -> Build the 3x3 linear system whose solution minimizes the quadric error
        //    associated with these two endpoints.
        // -> Use this system to solve for the optimal position, and store it in
        //    Edge_Record::optimal.
        // -> Also store the cost associated with collapsing this edge in
        //    Edge_Record::cost.
    }
    Halfedge_Mesh::EdgeRef edge;
    Vec3 optimal;
    float cost;
};

/* Comparison operator for Edge_Records so std::set will properly order them */
bool operator<(const Edge_Record& r1, const Edge_Record& r2) {
    if(r1.cost != r2.cost) {
        return r1.cost < r2.cost;
    }
    Halfedge_Mesh::EdgeRef e1 = r1.edge;
    Halfedge_Mesh::EdgeRef e2 = r2.edge;
    return &*e1 < &*e2;
}

/** Helper type for quadric simplification
 *
 * A PQueue is a minimum-priority queue that
 * allows elements to be both inserted and removed from the
 * queue.  Together, one can easily change the priority of
 * an item by removing it, and re-inserting the same item
 * but with a different priority.  A priority queue, for
 * those who don't remember or haven't seen it before, is a
 * data structure that always keeps track of the item with
 * the smallest priority or "score," even as new elements
 * are inserted and removed.  Priority queues are often an
 * essential component of greedy algorithms, where one wants
 * to iteratively operate on the current "best" element.
 *
 * PQueue is templated on the type T of the object
 * being queued.  For this reason, T must define a comparison
 * operator of the form
 *
 *    bool operator<( const T& t1, const T& t2 )
 *
 * which returns true if and only if t1 is considered to have a
 * lower priority than t2.
 *
 * Basic use of a PQueue might look
 * something like this:
 *
 *    // initialize an empty queue
 *    PQueue<myItemType> queue;
 *
 *    // add some items (which we assume have been created
 *    // elsewhere, each of which has its priority stored as
 *    // some kind of internal member variable)
 *    queue.insert( item1 );
 *    queue.insert( item2 );
 *    queue.insert( item3 );
 *
 *    // get the highest priority item currently in the queue
 *    myItemType highestPriorityItem = queue.top();
 *
 *    // remove the highest priority item, automatically
 *    // promoting the next-highest priority item to the top
 *    queue.pop();
 *
 *    myItemType nextHighestPriorityItem = queue.top();
 *
 *    // Etc.
 *
 *    // We can also remove an item, making sure it is no
 *    // longer in the queue (note that this item may already
 *    // have been removed, if it was the 1st or 2nd-highest
 *    // priority item!)
 *    queue.remove( item2 );
 *
 */
template<class T> struct PQueue {
    void insert(const T& item) {
        queue.insert(item);
    }
    void remove(const T& item) {
        if(queue.find(item) != queue.end()) {
            queue.erase(item);
        }
    }
    const T& top(void) const {
        return *(queue.begin());
    }
    void pop(void) {
        queue.erase(queue.begin());
    }
    size_t size() {
        return queue.size();
    }

    std::set<T> queue;
};

/*
    Mesh simplification. Note that this function returns success in a similar
    manner to the local operations, except with only a boolean value.
    (e.g. you may want to return false if you can't simplify the mesh any
    further without destroying it.)
*/
bool Halfedge_Mesh::simplify() {

    std::unordered_map<VertexRef, Mat4> vertex_quadrics;
    std::unordered_map<FaceRef, Mat4> face_quadrics;
    std::unordered_map<EdgeRef, Edge_Record> edge_records;
    PQueue<Edge_Record> edge_queue;

    // Compute initial quadrics for each face by simply writing the plane equation
    // for the face in homogeneous coordinates. These quadrics should be stored
    // in face_quadrics
    // -> Compute an initial quadric for each vertex as the sum of the quadrics
    //    associated with the incident faces, storing it in vertex_quadrics
    // -> Build a priority queue of edges according to their quadric error cost,
    //    i.e., by building an Edge_Record for each edge and sticking it in the
    //    queue. You may want to use the above PQueue<Edge_Record> for this.
    // -> Until we reach the target edge budget, collapse the best edge. Remember
    //    to remove from the queue any edge that touches the collapsing edge
    //    BEFORE it gets collapsed, and add back into the queue any edge touching
    //    the collapsed vertex AFTER it's been collapsed. Also remember to assign
    //    a quadric to the collapsed vertex, and to pop the collapsed edge off the
    //    top of the queue.

    // Note: if you erase elements in a local operation, they will not be actually deleted
    // until do_erase or validate are called. This is to facilitate checking
    // for dangling references to elements that will be erased.
    // The rest of the codebase will automatically call validate() after each op,
    // but here simply calling collapse_edge() will not erase the elements.
    // You should use collapse_edge_erase() instead for the desired behavior.

    return false;
}
