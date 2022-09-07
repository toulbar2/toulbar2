#ifndef __BST_h
#define __BST_h

/*
  Binary search Tree : AVL
*/

#include <iostream>

class Node {
public:
    Cost fitness;
    std::pair<int, int> neighbor;

    int height;
    Node * left;
    Node * right;

    Node(Cost _fitness, std::pair<int, int> & _neighbor) : fitness(_fitness), neighbor(_neighbor), height(1), left(NULL), right(NULL) {}

    Node(Node * _source) : fitness(_source->fitness), neighbor(_source->neighbor), height(1), left(NULL), right(NULL) {}

    void updateHeight() {
        if (left == NULL) {
            if (right == NULL)
                height = 1;
            else
                height = right->height + 1;
        } else
            if (right == NULL)
                height = left->height + 1;
            else
                height = std::max(left->height, right->height) + 1;
    }

    int balanceIndex() {
        if (left == NULL) {
            if (right == NULL)
                return 0;
            else
                return - right->height;
        } else
            if (right == NULL)
                return left->height ;
            else
                return left->height - right->height;
    }

    void print(std::ostream & _out) {
        _out << "{" << height << ", " << fitness << " " << neighbor.first << " " << neighbor.second << "}";
    }
};


class BST {
public:
    BST() : root(NULL) { }

    ~BST() {
        destroy(root);
        root = NULL;
    }

    void makeEmpty() {
        destroy(root);
        root = NULL;
    }

    void insert(Cost _fitness, std::pair<int, int> & _neighbor) {
        root = insert(root, _fitness, _neighbor);
    }

    void remove(Cost _fitness, std::pair<int, int> & _neighbor) {
        root = remove(root, _fitness, _neighbor);
    }

    void minimum(Cost & vmin, std::pair<int, int> & min) {
        Node * n = minimum(root);

        min.first  = n->neighbor.first;
        min.second = n->neighbor.second;
        vmin = n->fitness;
    }

    void findall(Cost _fitness, unsigned & _nb, std::vector< Node* > & _res) {
        _nb = 0;
        findall(root, _fitness, _nb, _res);
    }

    void print(std::ostream & _out) {
        prefixe(root, _out);
    }

    void copy(BST & _dest) {
        _dest.root = copy(root);
    }

    void copy(Node * _source, BST & _dest) {
        _dest.root = copy(_source);
    }

    bool checkHeight() {
        if (root == NULL)
            return true;
        else
            return (root->height == height(root));
    }

    bool checkBalance() {
        if (root == NULL)
            return true;
        else
            return checkBalance(root);
    }

    Node * root;

private:

    void destroy(Node * n) {
        if (n != NULL) {
            destroy(n->left);
            destroy(n->right);
            delete n;
        }
    }

    Node* copy(Node * n) {
        if (n == NULL) 
            return n;
        else {
            Node * cn = new Node(n);
            cn->left   = copy(n->left);
            cn->right  = copy(n->right);
            cn->height = n->height;

            return cn;
        }
    }

    int height(Node * n) {
        if (n == NULL)
            return 0;
        else {
            int l = height(n->left);
            int r = height(n->right);
            return std::max(l, r) + 1;
        }
    }

    bool checkBalance(Node * n) {
        if (n == NULL)
            return true;
        else {
            bool l = checkBalance(n->left);
            bool r = checkBalance(n->right);
            return l && r && (std::abs(n->balanceIndex()) < 2);
        }
    }

    Node * minimum(Node * n) {
        if (n == NULL)
            return NULL;
        else {
            if (n->left == NULL)
                return n;
            else
                return minimum(n->left);
        }
    }

    void findall(Node * n, Cost fit, unsigned & nb, std::vector< Node* > & res) {
        if (n != NULL) {
            if (fit == n->fitness) {
                res[nb] = n;
                nb++;
                findall(n->left, fit, nb, res);
                findall(n->right, fit, nb, res);                
            } else 
                if (fit < n->fitness) 
                    findall(n->left, fit, nb, res);
                else
                    findall(n->right, fit, nb, res);                
        }
    }

    Node* insert(Node * n, Cost _fitness, std::pair<int, int> & _neighbor) {
        if (n == NULL) {
            n = new Node(_fitness, _neighbor);
        } else {
            if (_fitness < n->fitness || 
                (_fitness == n->fitness && _neighbor.first < n->neighbor.first) || 
                (_fitness == n->fitness && _neighbor.first == n->neighbor.first && _neighbor.second < n->neighbor.second)) {
                n->left = insert(n->left, _fitness, _neighbor);
                n->updateHeight();
                n = balance(n);
            } else {
                n->right = insert(n->right, _fitness, _neighbor);
                n->updateHeight();
                n = balance(n);
            }
        }
        return n;
    }

    Node* remove(Node * n, Cost _fitness, std::pair<int, int> & _neighbor) {
        if (n != NULL) {
            if (n->fitness == _fitness && n->neighbor.first == _neighbor.first && n->neighbor.second == _neighbor.second) {
                // remove this node
                if (n->left == NULL) {
                    Node * tmp = n;

                    n = balance(n->right);
                    //n = balance(n);
                    delete tmp;
                } else {
                    std::pair<Node*, Node*> m = cutMax(n->left);
                    m.second->left  = m.first;

                    m.second->right = n->right;

                    m.second->updateHeight();
                    delete n;
                    n = m.second;
                }
            } else 
                if (_fitness < n->fitness || 
                    (_fitness == n->fitness && _neighbor.first < n->neighbor.first) || 
                    (_fitness == n->fitness && _neighbor.first == n->neighbor.first && _neighbor.second < n->neighbor.second)) {
                        n->left = remove(n->left, _fitness, _neighbor);
                        n->updateHeight();
                    } else {
                        n->right = remove(n->right, _fitness, _neighbor);
                        n->updateHeight();
                    }

            n = balance(n);
        }

        return n;
    }

    std::pair<Node*, Node*> cutMax(Node * n) {
        if (n->right == NULL) {
            return std::pair<Node*, Node*>(n->left, n);
        } else {
            std::pair<Node*, Node*> res = cutMax(n->right);

            n->right = res.first;

            n->updateHeight();
            
            n = balance(n);

            return std::pair<Node*, Node*>(n, res.second);            
        }
    }

    void prefixe(Node * n, std::ostream & _out) {
        _out << "(";
        if (n != NULL) {
            n->print(_out);
            prefixe(n->left, _out);
            prefixe(n->right, _out);
        }
        _out << ")";
    }

    int heightNode(Node * n) {
        return n ? n->height : 0;
    }

    Node* balance(Node * n) {
        if (n != NULL) {
            int b = n->balanceIndex();
            if (std::abs(b) > 1) {
                if (b == 2) {
                    if (heightNode(n->left->left) >= heightNode(n->left->right)) {
                        n = rotateRight(n);
                    } else {
                        n->left = rotateLeft(n->left);
                        n = rotateRight(n);
                    }
                } else {
                    if (heightNode(n->right->right) >= heightNode(n->right->left)) {
                        n = rotateLeft(n);
                    } else {
                        n->right = rotateRight(n->right);
                        n = rotateLeft(n);
                    }
                }
            }            
        }

        return n;
    }

    Node* rotateLeft(Node * p) {
        // right 
        Node * rt = p->right;

        // left of the right
        p->right = rt->left;
        p->updateHeight();

        rt->left = p;
        rt->updateHeight();

        return rt;
    }

    Node* rotateRight(Node * p) {
        // left
        Node * lt = p->left;

        p->left = lt->right;
        p->updateHeight();

        lt->right = p;
        lt->updateHeight();

        return lt;
    }
};

#endif
