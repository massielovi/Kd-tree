// Copyright

#ifndef SRC_KDTREE_HPP_
#define SRC_KDTREE_HPP_

#include <cmath>
#include <iostream>
#include <set>
#include <stdexcept>
#include <utility>
#include <vector>
#include <queue>
#include <unordered_map>
#include "Point.hpp"

using namespace std;

template <size_t N, typename ElemType>
class KDTree {
public:
    typedef pair<Point<N>, ElemType> value_type;
    KDTree();
    ~KDTree();
    KDTree(const KDTree& rhs);
    KDTree& operator=(const KDTree& rhs);
    size_t dimension() const;
    size_t size() const;
    bool empty() const;
    bool contains(const Point<N>& pt) const;
    void insert(const Point<N>& pt, const ElemType& value = ElemType());
    ElemType& operator[](const Point<N>& pt);
    ElemType& at(const Point<N>& pt);
    const ElemType& at(const Point<N>& pt) const;
    ElemType knn_value(const Point<N>& key, size_t k) const;
    vector<ElemType> knn_query(const Point<N>& key, size_t k) const;

private:
    size_t dimension_;
    size_t size_;

    struct KDNode
    {
        Point<N> pt_;
        KDNode* nodes[2];
        int level_;
        ElemType value_;

        KDNode(const Point<N>& pt, int level, const ElemType& value) :pt_(pt), level_(level), value_(value)
        {
            nodes[0] = nodes[1] = nullptr;
        }
    };

    KDNode* root;

    // TODO(me): finish the implementation of the rest of the KDTree class
    void nn_delete(KDNode* nnode);
    KDNode* copy(KDNode* root);
    KDNode* find_nn(const Point<N>& pt, KDNode* node_a) const;
    //void nn_Search(const typename KDNode* node_c, const Point<N>& key, priority_queue<ElemType>& queque, size_t k) const ;

};

// TODO(me):the KDTree class

template <size_t N, typename ElemType>
KDTree<N, ElemType>::KDTree() {
    // Inicializamos N=3
    this->root = nullptr;
    size_ = 0;
}

template <size_t N, typename ElemType>
KDTree<N, ElemType>::~KDTree() {
    /*delete root->nodes[0];delete root->nodes[1];
    delete this->root;
    */
    nn_delete(this->root);
}

template <size_t N, typename ElemType>
KDTree<N, ElemType>::KDTree(const KDTree& rhs) {
    this->root = copy(rhs.root);
    this->size_ = rhs.size_;
}

template <size_t N, typename ElemType>
KDTree<N, ElemType>& KDTree<N, ElemType>::operator=(const KDTree& rhs) {
    if (this != &rhs)
    {
        nn_delete(this->root);
        this->root = copy(rhs.root);
        this->size_ = rhs.size_;
    }
    return *this;
}

template <size_t N, typename ElemType>
size_t KDTree<N, ElemType>::dimension() const {
    return N;
}

template <size_t N, typename ElemType>
size_t KDTree<N, ElemType>::size() const {
    return size_;
}

template <size_t N, typename ElemType>
bool KDTree<N, ElemType>::empty() const {
    return size_ == 0;
}

template <size_t N, typename ElemType>
bool KDTree<N, ElemType>::contains(const Point<N>& pt) const {
    KDNode* node = find_nn(pt, this->root);

    if (node != nullptr && node->pt_ == pt)return true;

    return false;
}

template <size_t N, typename ElemType>
void KDTree<N, ElemType>::insert(const Point<N>& pt, const ElemType& value) {

    KDNode* node_b = find_nn(pt, this->root);
    if (node_b == nullptr)
    {
        this->root = new KDNode(pt, 0, value);
        this->size_ = 1;
    }
    else
    {
        if (node_b->pt_ == pt) node_b->value_ = value;
        else
        {
            int Level_b = node_b->level_;
            KDNode* newNode = new KDNode(pt, Level_b + 1, value);

            if (pt[Level_b % N] < node_b->pt_[Level_b % N]) node_b->nodes[0] = newNode;
            else node_b->nodes[1] = newNode;

            ++size_;
        }
    }

}

template <size_t N, typename ElemType>
ElemType& KDTree<N, ElemType>::operator[](const Point<N>& pt) {

    KDNode* node = find_nn(pt, this->root);

    if (node != nullptr && node->pt_ == pt) return node->value_;
    else
    {
        insert(pt);

        if (node == NULL) return this->root->value_;
        else
        {
            if (node->nodes[0] != NULL && node->nodes[0]->pt_ == pt)return node->nodes[0]->value_;
            else return node->nodes[1]->value_;
        }
    }
}

template <size_t N, typename ElemType>
ElemType& KDTree<N, ElemType>::at(const Point<N>& pt) {
    const KDTree<N, ElemType>& constThis = *this;
    return const_cast<ElemType&>(constThis.at(pt));
}

template <size_t N, typename ElemType>
const ElemType& KDTree<N, ElemType>::at(const Point<N>& pt) const {
    KDNode* node = find_nn(pt, this->root);

    if (node == nullptr || node->pt_ != pt) throw out_of_range("No existe Point");
    else return node->value_;

}
/*
template <size_t N, typename ElemType>
ElemType KDTree<N, ElemType>::knn_value(const Point<N> &key, size_t k) const {
  // TODO(me): Fill this in.
  priority_queue<ElemType> Qqueque(k);

  if (empty()) return ElemType(); is empty


    nn_Search(root, key, Qqueque,k);


    unordered_map<ElemType, int> counter;
    int tmp = Qqueque.begin()->second;
    Qqueque.erase(Qqueque.begin());

    while (!Qqueque.empty())
      {
        ++counter[tmp];
    }
    ElemType new_element;

    int cnt = -1;
    for (const auto &p : counter)
    {
        if (p.second > cnt)
        {
            new_element = p.first;
            cnt = p.second;
        }
    }
  return new_element;
}
*/

template <size_t N, typename ElemType>
vector<ElemType> KDTree<N, ElemType>::knn_query(const Point<N>& key,
    size_t k) const {
    // TODO(me): Fill this in.
    vector<ElemType> values;
    return values;
}

// TODO(me): finish the implementation of the rest of the KDTree class


template <size_t N, typename ElemType>
void KDTree<N, ElemType>::nn_delete(typename KDTree<N, ElemType>::KDNode* nnode)
{
    if (nnode == nullptr) return;

    nn_delete(nnode->nodes[0]); nn_delete(nnode->nodes[1]);

    --size_;
    delete nnode;
}

template <size_t N, typename ElemType>
typename KDTree<N, ElemType>::KDNode* KDTree<N, ElemType>::copy(typename KDTree<N, ElemType>::KDNode* root_)
{
    if (root_ == nullptr) return NULL;

    KDNode* copy_root = new KDNode(*root_);

    copy_root->nodes[0] = copy(copy_root->nodes[0]);
    copy_root->nodes[1] = copy(copy_root->nodes[1]);

    return copy_root;
}

template <size_t N, typename ElemType>
typename KDTree<N, ElemType>::KDNode* KDTree<N, ElemType>::find_nn(const Point<N>& pt, typename KDTree<N, ElemType>::KDNode* node_a) const
{
    if (node_a == NULL || node_a->pt_ == pt) return node_a;

    const Point<N>& Point_a = node_a->pt_;
    int Level_a = node_a->level_;

    if (pt[Level_a % N] < Point_a[Level_a % N])
    {
        if (node_a->nodes[0] == NULL)return node_a;
        else return find_nn(pt, node_a->nodes[0]);
    }
    else
    {
        if (node_a->nodes[1] == NULL)return node_a;
        else return find_nn(pt, node_a->nodes[1]);
    }
}

/*
template <std::size_t N, typename ElemType>
void KDTree<N, ElemType>::nn_Search(const typename KDTree<N, ElemType>::KDNode* node_c, const Point<N>& key, priority_queue<ElemType>& queque, size_t k) const
{
    if (node_c == nullptr) return;
    const Point<N>& Point_c = node_c->point;
    queque.push(pair(Distance(Point_c, key),node_c->value));

    if (queque.size() > k) {
        typename multimap<double, T>::iterator last = queque.end();
        --last;
        queque.erase(last);
    }

    int Level_c = node_c->level;
    bool isLeftTree;
    if (key[Level_c%N] < Point_c[Level_c%N])
    {
        nn_Search(node_c->left, key, queque,k);
        isLeftTree = true;
    }
    else
    {
        nn_Search(node_c->right, key, queque,k);
        isLeftTree = false;
    }

    double worst;
    if(empty()) worst=numeric_limits<double>::infinity() ;
    else worst=queque.begin()->first;

    if (queque.size() < k || fabs(key[Level_c%N] - Point_c[Level_c%N]) < worst)
    {
        if (isLeftTree) nn_Search(node_c->right, key, queque,k);
        else nn_Search(node_c->left, key, queque,k);
    }
}
*/
#endif  // SRC_KDTREE_HPP_
