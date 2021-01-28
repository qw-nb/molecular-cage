#include <iostream>
#include <vector>
using namespace std;
/* the type of rotation group */
enum group{ C=0, D=1, T=2, O=3, I=4 };
/* the type of element */
enum element{ Face=0, Window=1, Edge=2, Vertex=3 };
/* the face element */
struct face
{
    int degree; // the degree of the face element
    int num;    // the number of the face
};
/* the vertex of polyhedron */
struct vertex
{
    int degree; // normally only 3 and 4 is permitted
    int num;    // the number of the vertex
};
/* the axis of polyhedron */
struct axis
{
    int degree; // the order of axis
    int num;    // the number of axis
    element e_1, e_2;   // the element passed by axis
};
/* the connection between nodes */
struct connection
{
    int layer_1, layer_2; // the order of node on the chains
    int chain_number;  // if the int is 0 if the connection is on the same chain, else on different ones
};
struct node
{
    int deep;   // the deepth of nodes in tree
    connection data;   // store the position of data in dataset 
    node* parent_node;
    vector<node*> son_node;
};
void get_allconn(int chain_len, vector<connection>& all_iconn, vector<connection>& all_oconn);
void traverse(node* nd, int deep, vector<connection*>& result, connection temp[]);
vector<vector<int>> get_allcycle(node* nd, int chain_len,const vector<connection>& all_oconn,const vector<connection>& all_iconn);
bool if_matchc(vector<vector<int>> test, face f, face window, bool iff1, bool iff2, int chain_len);
bool if_matchv(vector<connection> test, int chain_len, int vn[], vertex v[], connection iff1, connection iff2);
bool if_macthp(const vector<int>& a,const vector<int>& b, int chain_len);
void node_delete(node* nd);
void clean_data(vector<connection>& dataset, connection data);
void make_tree(node* root, int max_deep, vector<connection>& dataset, int i, connection iff1,connection iff2,  
                vector<connection> iconn, vector<connection> oconn,face f, face w, int vn[], vertex v[],int chain_len);
void grow_allleaf(node* root, int max_deep, vector<connection> data,connection iff1,connection iff2,  
            vector<connection> iconn,face f, face w, int vn[], vertex v[],int chain_len);
/* the class of polyhedron */
class polyhedron{
    private:
    face f; // the face element of polyhedron
    face w; // the window element of polyhedron
    vertex* v;  // the vertex of polyhedron, the len is 2
    /* the axises of the polyhedron 
        C group with only one main axis
        D group the first is C2 axis, the second is main axis
        T group the first is C2 axis, the second is C3 axis
        O group the first is C2 axis, the second is C3 axis, the last is C4 axis
        I group the first is C2 axis, the second is C3 axis, the last is C5 axis */
    vector<axis> a; 
    group type; // the group of the polyhedron
    int chain_len;  // the length of chain
    int conn_num;   // the number of connection on one chain
    polyhedron(){}; // forbide using of auto struct function
    public:
    polyhedron(face f_, face w_, vector<axis> a_, group t):
    f(f_), w(w_), a(a_), type(t){
        v=new vertex[2];
        v[0].degree=3;
        v[1].degree=4;
        int edge=f.num+w.num+(f.num*f.degree)/2-2; // the number of edges
        int all_vn=f.num*f.degree/2;    // the total number of vertexs
        int v4n=edge+edge-all_vn*3;   // the number of vertex of 4 degree
        int v3n=all_vn-v4n; // the number of vertex of 3 degree
        v[0].num=v3n; v[1].num=v4n;
        int order=0;    // the order of group
        for (int i=0; i<a.size(); i++){
            switch (a[i].e_1)
            {
            case Face:
                f.num=f.num-a[i].num;
                break;
            case Window:
                w.num=w.num-a[i].num;
                break;
            case Edge:
                edge=edge-a[i].num;
                break;
            case Vertex:
                v[1].num=v[1].num-a[i].num;
                break;
            default:
                break;
            }
            switch (a[i].e_2)
            {
            case Face:
                f.num=f.num-a[i].num;
                break;
            case Window:
                w.num=w.num-a[i].num;
                break;
            case Edge:
                edge=edge-a[i].num;
                break;
            case Vertex:
                v[1].num=v[1].num-a[i].num;
                break;
            default:
                break;
            }
        }
        /* determine the order of the polyhedron's group */ 
        switch (type)
        {
        case C:
            order=a[0].degree;
            break;
        case D:
            order=2*a[0].degree;
            break;
        case T:
            order=12;
            break;
        case O:
            order=24;
            break;
        case I:
            order=60;
            break;
        default:
            break;
        }
        chain_len=(v[0].num+v[1].num)/order;
        conn_num=edge/order-chain_len+1;
    };
    polyhedron(const polyhedron& orther):
    f(orther.f), w(orther.w), v(orther.v), a(orther.a),
    chain_len(orther.chain_len), conn_num(orther.conn_num),type(orther.type)
    {};
    polyhedron& operator=(const polyhedron& orther){
        f=orther.f;
        w=orther.w;
        a=orther.a;
        v=orther.v;
        chain_len=orther.chain_len;
        conn_num=orther.conn_num;
        type=orther.type;
        return *this;
    }
    int get_chainlen(){
        return chain_len;
    }
    int get_connnum(){
        return conn_num;
    }
    vertex* get_v(){
        return v;
    }
};
class tree{
    friend void traverse(node* nd, int deep, vector<int*> result, int temp[]);
    private:
    int deep;   // the maxium of deepth for tree
    node* root;
    vector<connection> idata_set;   // store in-connection datas
    vector<connection> odata_set;   // store out-connection datas
    tree(){}
    tree(const tree& orther){}
    public:
    tree(int d,vector<connection> ids, vector<connection> ods):
    deep(d), idata_set(ids),odata_set(ods){
        root=new node;
        root->deep=0;
        root->data={-1,-1,0};
        root->parent_node=nullptr;
    }
    void clean_dataset(face f, face w){
        int fd=f.degree, wd=w.degree;
        int i=0;
        while(i<idata_set.size()){
            if (idata_set[i].chain_number){
                int rn=abs(idata_set[i].layer_2-idata_set[i].layer_1)+1;
                if (rn!=fd && rn!=wd){
                    idata_set.erase(idata_set.begin()+i);
                }
                else i++;
            }
            else i++;
        }
    }
    void grow_Ctree(const vector<axis>& a, int chain_len, face f, face w, vertex v[]){
        int vn[chain_len];  // judge if the connection is over degree
        node* pointer=root;
        connection iff1{-1,-1,0},iff2{-2,-2,0};
        vector<connection> all_iconn, all_oconn;
        for (int i=0; i<chain_len; i++){
            vn[i]=2;
        }
        vn[0]--; vn[chain_len-1]--;
        if (a[0].e_1==Vertex){
            vn[0]++; vn[1]++;
            iff1={0,f.degree-3,false};
            all_oconn.push_back(iff1);
            w.num--; f.num--;
        }
        else {
            iff1={0,0,false};
            if (a[0].e_1==Face){
            }
        }
        if (a[0].e_2==Vertex){
            vn[0]++; vn[1]++;
            iff2={chain_len-1,chain_len-f.degree+2,false};
            all_oconn.push_back(iff2);
            w.num--; f.num--;
        }
        else {
            iff2={chain_len-1,chain_len-1,false};
            if (a[0].e_1==Face){
            }
        }
        clean_data(odata_set,iff1); clean_data(odata_set,iff2);
        pointer->data=iff1;
        node temp{1,iff2,root};
        node* u=&temp;
        u->parent_node=pointer;
        pointer->son_node.push_back(u); pointer=pointer->son_node[0];
        vector<connection> iconn;
        vector<connection> oconn;
        for (int i=2; i<=deep; i++){
            make_tree(pointer,i,idata_set,0,iff1,iff2,iconn,oconn,f,w,&vn[0],&v[0],chain_len);
            cout <<"sd" <<root->son_node[0]->son_node.size() << "  "<< idata_set.size()<< endl;
            cout << i << endl;
        }
        grow_allleaf(pointer,deep,odata_set,iff1,iff2,iconn,f,w,vn,v,chain_len);
        vector<connection*> result; connection t[deep];
        traverse(root,deep-1,result,t);

    }
    void get_root(){
        cout << "hellp";
    }
};
void traverse(node* nd, int deep, vector<connection*>& result, connection temp[]){
    if (deep>0){
        temp[deep-1]=nd->data;
        for (int i=0; i<nd->son_node.size(); i++){
            return traverse(nd->son_node[i], deep-1, result, temp);
        }
    }
    else{
        cout << "sd";
        temp[deep]=nd->data;
        result.push_back(&temp[0]);
        cout << result.size() << "ese" << endl;
        return;
    }
}
vector<vector<int>> get_allcycle(node* nd, int chain_len,const vector<connection>& all_oconn,
                    const vector<connection>& all_iconn, connection axis1, connection axis2){
    vector<vector<int>> result;
    for (int i=0; i<all_iconn.size(); i++){
        vector<int> temp;
        for (int j=all_iconn[i].layer_1; j<=all_iconn[i].layer_2; j++){
            temp.push_back(j);
        }
        result.push_back(temp);
    }
    int i=0; vector<int> temp; bool if_path=false;
    if (all_oconn.size()<=0) return result;
    for (int j=0; j<all_iconn.size(); j++){
        if (all_iconn[j].layer_1>=axis1.layer_1 && all_iconn[j].layer_2<=all_oconn[0].layer_1){
            temp.push_back(axis1.layer_1);
            for (int k=axis1.layer_1; k<=all_iconn[j].layer_1; k++){
                temp.push_back(k);
            }
            for (int k=all_iconn[j].layer_2; k<=all_oconn[i].layer_1; k++){
                temp.push_back(k);
            }
            for (int k=all_oconn[i].layer_2+chain_len; k>=chain_len+axis1.layer_1; k--){
                temp.push_back(k);
            }
            if_path=true;
            break;
        }
        else if (all_iconn[j].layer_1>=axis1.layer_2 && all_iconn[j].layer_2<=all_oconn[0].layer_2){
            for (int k=axis1.layer_1; k<=all_oconn[i].layer_1; k++){
                temp.push_back(k);
            }
            for (int k=all_oconn[i].layer_2+chain_len; k>=all_iconn[j].layer_2+chain_len; k--){
                temp.push_back(k);
            }
            for (int k=all_oconn[j].layer_1+chain_len; k>=axis1.layer_2+chain_len; k--){
                temp.push_back(k);
            }
            if_path=true;
            break;
        }
    }
    if (!if_path){
        for (int j=axis1.layer_1; j<=all_oconn[i].layer_1;j++){
            temp.push_back(j);
        }
        for (int j=all_oconn[i].layer_2+chain_len; j>=axis1.layer_2+chain_len; j--){
            temp.push_back(j);
        }
    }
    result.push_back(temp);
    for (i=0; i<all_oconn.size(); i++){
        temp.clear(); if_path=false;
        if (i==all_oconn.size()-1){
            for (int j=0; j<all_iconn.size(); j++){
                if (all_iconn[j].layer_1>=all_oconn[i].layer_1 && all_oconn[j].layer_2<=axis2.layer_1){
                    for (int k=all_oconn[i].layer_1; k<=all_iconn[j].layer_1; k++){
                        temp.push_back(k);
                    }
                    for (int k=all_iconn[j].layer_2; k<=axis2.layer_1; k++){
                        temp.push_back(k);
                    }
                    for (int k=axis2.layer_2+chain_len; k>=all_oconn[i].layer_2+chain_len; k--){
                        temp.push_back(k);
                    }
                    if_path=true;
                    break;
                }
                else if (all_iconn[j].layer_1>=all_oconn[i].layer_2 && all_oconn[j].layer_2<=axis2.layer_2){
                    for (int k=all_oconn[i].layer_1; k<=axis2.layer_1; k++){
                        temp.push_back(k);
                    }
                    for (int k=axis2.layer_2+chain_len; k>=all_iconn[j].layer_2+chain_len; k--){
                        temp.push_back(k);
                    }
                    for (int k=all_iconn[j].layer_1+chain_len; k>=all_oconn[i].layer_2+chain_len; k--){
                        temp.push_back(k);
                    }
                    if_path=true;
                    break;
                }
            }
            if (!if_path){
                for (int j=all_oconn[i].layer_1; j<=axis2.layer_1; j++){
                    temp.push_back(j);
                }
                for (int j=axis2.layer_2+chain_len; j>=all_oconn[i].layer_2+chain_len; j--){
                    temp.push_back(j);
                }
            }
        }
        else {
            for (int j=0; j<all_iconn.size(); j++){
                if (all_iconn[j].layer_1>=all_oconn[i].layer_1 && all_iconn[j].layer_2<=all_oconn[i+1].layer_1){
                    for (int k=all_oconn[i].layer_1; k<=all_iconn[j].layer_1; k++){
                        temp.push_back(k);
                    }
                    for (int k=all_iconn[j].layer_2; k<=all_oconn[i+1].layer_1; k++){
                        temp.push_back(k);
                    }
                    for (int k=all_oconn[i+1].layer_2+chain_len; k>=all_oconn[i].layer_2+chain_len; k--){
                        temp.push_back(k);
                    }
                    if_path=true;
                    break;
                }
                else if (all_iconn[j].layer_1>=all_oconn[i].layer_2 && all_iconn[j].layer_2<=all_oconn[i+1].layer_2){
                    for (int k=all_oconn[i].layer_1; k<=all_oconn[i+1].layer_1; k++){
                        temp.push_back(k);
                    }
                    for (int k=all_oconn[i+1].layer_2+chain_len; k>=all_iconn[j].layer_2; k--){
                        temp.push_back(k);
                    }
                    for (int k=all_iconn[j].layer_1; k>=all_oconn[i].layer_2; k--){
                        temp.push_back(k);
                    }
                    if_path=true;
                    break;
                }
            }
            if (!if_path){
                for (int j=all_oconn[i].layer_1; j<=all_oconn[i+1].layer_1; j++){
                    temp.push_back(j);
                }
                for (int j=all_oconn[i+1].layer_2+chain_len; j>=all_oconn[i].layer_2+chain_len; j--){
                    temp.push_back(j);
                }
            }
        }
        result.push_back(temp);
    }
    return result;
}
void get_allconn(int chain_len, vector<connection>& all_iconn, vector<connection>& all_oconn){
    for (int i=0; i<chain_len; i++){
        for (int j=0; j<i; j++){
            connection temp={j,i,true};
            all_iconn.push_back(temp); // make the connection in one chain
        }
    }
    for (int i=0; i<chain_len; i++){
        for (int j=0; j<chain_len; j++){
            connection temp={i,j,false};
            all_oconn.push_back(temp); // make the connection between the chains
        }
    }
    return;
}
/*
    to judge if the combination of connections fits the vertexs' need
    vn is the max degree the node can be, chain_len is this dynamic list's length
*/
bool if_matchv(vector<connection> test, int chain_len, int vn[], vertex v[], connection iff1, connection iff2){
    int vn_[chain_len];
    int v3n=0, v4n=0;
    for (int i=0; i<chain_len; i++){
        vn_[i]=vn[i];
    }
    vn_[iff1.layer_1]++; vn_[iff1.layer_2]++;
    vn_[iff2.layer_1]++; vn_[iff2.layer_2]++;
    for (int i=0; i<test.size(); i++){
        vn_[test[i].layer_1]++;
        vn_[test[i].layer_2]++;
    }
    for (int i=0; i<chain_len; i++){
        if (vn_[i]>4) return false;
        else if (vn_[i]==4) v4n++;
        else if (vn_[i]==3) v3n++;
    }
    if (v4n>v[0].num || v3n>v[1].num) {
        return false;
    }
    return true;
}
bool if_matchc(vector<vector<int>> test, face f, face window, connection iff1, connection iff2, int chain_len){
    vector<vector<int>> fcycle;
    vector<vector<int>> wcycle;
    for (int i=0; i<test.size()-1; i++){
        if (test[i].size()==f.degree){
            fcycle.push_back(test[i]);
        }
        else if (test[i].size()==window.degree){
            if (test[i][0]==iff1.layer_1 && test[i].back()==iff1.layer_2+chain_len) {
                return false;
            }
            else if (test[i][0]==iff2.layer_1 && test[i].back()==iff2.layer_2+chain_len) {
                return false;
            }
            for (int j=0; j<wcycle.size(); j++){
                if (if_macthp(test[i],wcycle[j],chain_len)) {
                    return false;
                }
            }
            wcycle.push_back(test[i]);
        }
        else return false;
    }
    if (f.degree!=window.degree){
        if (fcycle.size()<=f.num && wcycle.size()<=window.num) return true;
    }
    else {
        if (fcycle.size() + wcycle.size()<=window.num+f.num) return true;
    }
    return false;
}
bool if_macthp(const vector<int>& a,const vector<int>& b, int chain_len){
    bool result=false;
    for (int i=0; i<a.size(); i++){
        for (int j=0; j<b.size(); j++){
            if (a[i]==b[j] && a[i+1]==b[j+1]){
                result=true; break;
            }
            else if (a[i]+chain_len==b[j] && a[i+1]+chain_len==b[j+1]){
                result=true; break;
            }
            else if (a[i+1]==b[j+1]+chain_len && a[i]==b[j]+chain_len){
                result=true; break;
            }
        }
   }
   return result;
}
void clean_data(vector<connection>& dataset, connection data){
    for (int i=0; i<dataset.size(); i++){
        if (dataset[i].layer_1==data.layer_1 && dataset[i].layer_2==data.layer_2){
            dataset.erase(dataset.begin()+i);
            return;
        }
    }
}
void node_insert(node* pnd, node snd){
    node* pointer=&snd;
    pnd->son_node.push_back(pointer);
    return;
}
void node_delete(node* nd){
    node* pointer=nd->parent_node;
    if (pointer!=nullptr){
        for (int i=0; i<pointer->son_node.size(); i++){
            if (pointer->son_node[i]->data.layer_1==nd->data.layer_1 &&
             pointer->son_node[i]->data.layer_2==nd->data.layer_2){
                pointer->son_node.erase(pointer->son_node.begin()+i);
                return;
            }
        }
    }
}
void make_tree(node* root, int max_deep, vector<connection>& dataset, int i, connection iff1,connection iff2,  
                vector<connection> iconn, vector<connection> oconn,face f, face w, int vn[], vertex v[],int chain_len){ 
    if (root->deep<max_deep){
        int s=dataset.size()-max_deep+root->deep+1;
        node* u=new node[s];
        for (int j=i; j<=dataset.size()-max_deep+root->deep; j++){
            node* pointer=root;
            vector<connection> iconn_=iconn, oconn_=oconn;
            if (dataset[j].chain_number) iconn_.push_back(dataset[j]);
            else {oconn_.push_back(dataset[j]);}
            vector<vector<int>> cycle=get_allcycle(root,chain_len,oconn_,iconn_,iff1,iff2);
            vector<connection> vtest=iconn_;
            for (int k=0; k<oconn_.size(); k++){
                vtest.push_back(oconn_[k]);
            }
            if ((!if_matchc(cycle,f,w,iff1,iff2,chain_len))||(!if_matchv(vtest,chain_len,vn,v,iff1,iff2))){

                continue;
            }
            node temp{root->deep+1,dataset[j],root};
            u[j]=temp;
            root->son_node.push_back(&u[j]);
            make_tree(&u[j],max_deep,dataset,j+1,iff1,iff2,iconn_,oconn_,f,w,vn,v,chain_len);
        }
        
    }
    else {
        return;
    }
}
void grow_allleaf(node* root, int max_deep, vector<connection> data,connection iff1,connection iff2,  
            vector<connection> iconn,face f, face w, int vn[], vertex v[],int chain_len){
    if (root->son_node.size()==0){
        vector<connection> oconn;
        make_tree(root,max_deep,data,0,iff1,iff2,iconn,oconn,f,w,vn,v,chain_len);
    }
    else{
        for (int i=0; i<root->son_node.size(); i++){
            vector<connection> iconn_=iconn;
            iconn_.push_back(root->son_node[i]->data);
            grow_allleaf(root->son_node[i],max_deep,data,iff1,iff2,iconn_,f,w,vn,v,chain_len);
        }
    }
}
