#ifndef TSPSOLVER_TOUR_TREE_HPP
#define TSPSOLVER_TOUR_TREE_HPP

#include <utility>
#include <vector>
#include <stack>
#include <unordered_set>
#include <iostream>

using namespace std;

struct SplayNode {
    SplayNode* left = nullptr;
    SplayNode* right = nullptr;
    SplayNode* parent = nullptr;
    int value = 0;
    bool reversed = false;
};

struct SplayTree {
    SplayNode* root = nullptr;
    vector<SplayNode*> value_to_node;
    int size = 0;

    explicit SplayTree(const vector<int>& tour_sequence) {
        value_to_node = vector<SplayNode*>(tour_sequence.size());
        initialize_with_tour_seq(tour_sequence);
    };

    explicit SplayTree(SplayTree* splayTree) {
        vector<int>* sequence = splayTree->get_sequence();
        initialize_with_tour_seq(*sequence);
        delete sequence;
    };

    ~SplayTree() {
        stack<SplayNode*> stack_of_nodes;
        stack_of_nodes.push(root);
        while (!stack_of_nodes.empty()) {
            SplayNode* current_node = stack_of_nodes.top();
            if (current_node->left) {
                stack_of_nodes.push(current_node->left);
            }
            if (current_node->right) {
                stack_of_nodes.push(current_node->right);
            }
            stack_of_nodes.pop();
            delete current_node;
        }
    }

    void initialize_with_tour_seq(const vector<int>& tour_sequence) {
        for (auto value : tour_sequence) {
            this->insert(value);
        }
    }

    static void rightRotate(SplayNode *ptr) {
        SplayNode *left=ptr->left;
        SplayNode *left_right=left->right;
        SplayNode *parent=ptr->parent;
        if(parent) {
            if(parent->right == ptr) {
                parent->right=left;
            } else {
                parent->left = left;
            }
        }
        if(left_right) {
            left_right->parent = ptr;
        }
        left->parent=parent;
        left->right=ptr;

        ptr->parent=left;
        ptr->left=left_right;
    }

    static void leftRotate(SplayNode *ptr) {
        SplayNode *right=ptr->right;
        SplayNode *right_left=right->left;
        SplayNode *parent=ptr->parent;
        if(parent) {
            if(parent->right == ptr) {
                parent->right=right;
            } else {
                parent->left=right;
            }
        }
        if(right_left) {
            right_left->parent = ptr;
        }
        right->parent=parent;
        right->left=ptr;

        ptr->parent=right;
        ptr->right=right_left;
    }

    static void cleanUpNodeReversal(SplayNode *ptr, SplayNode *p, SplayNode *pp) {
        passDownReversal(pp);
        passDownReversal(p);
        passDownReversal(ptr);
    }

    static void passDownReversal(SplayNode *ptr) {
        if (ptr != nullptr && ptr->reversed) {
            ptr->reversed=false;
            swap(ptr->left, ptr->right);

            if(ptr->left != nullptr) {
                ptr->left->reversed = !ptr->left->reversed;
            }
            if(ptr->right != nullptr) {
                ptr->right->reversed = !ptr->right->reversed;
            }
        }
    }

    void Splay(SplayNode *key) {
        while(true) {
            SplayNode *p=key->parent;
            if(!p) {
                break;
            }
            SplayNode *pp=p->parent;

            cleanUpNodeReversal(key, p, pp);

            if (!pp) {
                if(p->left == key) {
                    rightRotate(p);
                } else {
                    leftRotate(p);
                }
                break;
            }
            if (pp->left==p) {
                if(p->left == key) {
                    rightRotate(pp);
                    rightRotate(p);
                } else {
                    leftRotate(p);
                    rightRotate(pp);
                }
            } else {
                if(p->left == key) {
                    rightRotate(p);
                    leftRotate(pp);
                } else {
                    leftRotate(pp);
                    leftRotate(p);
                }
            }
        }
        root=key;
    }
    
    SplayNode* insert(int value) {
        auto* node = new SplayNode();
        node->value = value;
        value_to_node[value] = node;

        if(!root) {
            root = node;
            size++;
            return node;
        }
        SplayNode* ptr = root;
        while (ptr->right != nullptr) {
            ptr = ptr->right;
        }
        ptr->right = node;
        node->parent = ptr;
        size++;
        return node;
    }

    [[nodiscard]] vector<int>* get_sequence() const {
        auto* sequence = new vector<int>(size);
        int current_index_number = 0;
        stack<pair<SplayNode*, bool>> stack_of_nodes;
        SplayNode* current_ptr = root;

        bool reversed = false;

        while (current_ptr != nullptr || !stack_of_nodes.empty()) {
            while(current_ptr != nullptr) {
                if (current_ptr->reversed) {
                    reversed = !reversed;
                }
                stack_of_nodes.push(make_pair(current_ptr, reversed));
                if (reversed) {
                    current_ptr = current_ptr->right;
                } else {
                    current_ptr = current_ptr->left;
                }
            }

            auto stacked_item = stack_of_nodes.top();
            current_ptr = stacked_item.first;
            reversed = stacked_item.second;
            stack_of_nodes.pop();

            (*sequence)[current_index_number] = current_ptr->value;
            current_index_number++;

            if (reversed) {
                current_ptr = current_ptr->left;
            } else {
                current_ptr = current_ptr->right;
            }
        }

        return sequence;
    }

    SplayNode* find(int value) {
        if (value < 0 || value > value_to_node.size()) {
            return nullptr;
        }

        SplayNode* node = value_to_node[value];
        Splay(node);
        return node;
    }

    void flip(int b, int d) {
        if (b == d) {
            return;
        }
        find(d);
        find(b);

        if (root->left != nullptr && root->left->value == d) {
            swapNodesLeft();
            if (root->left->right != nullptr) {
                root->left->right->reversed = !root->left->right->reversed;
            }
            root->reversed = !root->reversed;
        } else if (root->left != nullptr && root->left->left != nullptr && root->left->left->value == d) {
            rightRotate(root->left);
            swapNodesLeft();
            if (root->left->right != nullptr) {
                root->left->right->reversed = !root->left->right->reversed;
            }
            root->reversed = !root->reversed;
        } else if (root->right != nullptr && root->right->value == d) {
            if (root->right->left != nullptr) {
                root->right->left->reversed = !root->right->left->reversed;
            }
        } else if (root->right != nullptr && root->right->right !=nullptr && root->right->right->value == d) {
            leftRotate(root->right);
            if (root->right->left != nullptr) {
                root->right->left->reversed = !root->right->left->reversed;
            }
        } else {
            cout << "WEIRD CASE";
        }
    }

    void swapNodesLeft() {
        SplayNode* bNode = root->left;
        SplayNode* bNodeLeft = bNode->left;
        SplayNode* bNodeRight = bNode->right;
        SplayNode* parentRight = root->right;
        bNode->parent = nullptr;
        bNode->left = root;
        bNode->right = parentRight;
        root->left = bNodeLeft;
        root->right = bNodeRight;
        root->parent = bNode;
        if (bNodeLeft != nullptr) {
            bNodeLeft->parent = root;
        }
        if (bNodeRight != nullptr) {
            bNodeRight->parent = root;
        }
        if (parentRight != nullptr) {
            parentRight->parent = bNode;
        }
        root = bNode;
    }

    int nextValue(int value) {
        SplayNode* valueNode = value_to_node[value];
        Splay(valueNode);

        if (root == nullptr) {
            return 0;
        }

        if (root->left == nullptr && root->right == nullptr) {
            return root->value;
        }

        SplayNode* ptr;
        if (root->right == nullptr) {
            ptr = root->left;
        } else {
            ptr = root->right;
        }
        bool reversed = false;
        while(true) {
            if (ptr->reversed) {
                reversed = !reversed;
            }
            if (reversed) {
                if (ptr->right != nullptr) {
                    ptr = ptr->right;
                } else {
                    break;
                }
            } else {
                if (ptr->left != nullptr) {
                    ptr = ptr->left;
                } else {
                    break;
                }
            }
        }
        return ptr->value;
    }

    int previousValue(int value) {
        SplayNode* valueNode = value_to_node[value];
        Splay(valueNode);

        if (root == nullptr) {
            return 0;
        }

        if (root->left == nullptr && root->right == nullptr) {
            return root->value;
        }

        SplayNode* ptr;
        if (root->left == nullptr) {
            ptr = root->right;
        } else {
            ptr = root->left;
        }
        bool reversed = false;
        while(true) {
            if (ptr->reversed) {
                reversed = !reversed;
            }
            if (reversed) {
                if (ptr->left != nullptr) {
                    ptr = ptr->left;
                } else {
                    break;
                }
            } else {
                if (ptr->right != nullptr) {
                    ptr = ptr->right;
                } else {
                    break;
                }
            }
        }
        return ptr->value;
    }

    void swapNodesRight() {
        SplayNode* bNode = root->right;
        SplayNode* bNodeLeft = bNode->left;
        SplayNode* bNodeRight = bNode->right;
        SplayNode* parentLeft = root->left;
        bNode->parent = nullptr;
        bNode->right = root;
        bNode->left = parentLeft;
        root->parent = bNode;
        root->left = bNodeLeft;
        root->right = bNodeRight;
        if (bNodeLeft != nullptr) {
            bNodeLeft->parent = root;
        }
        if (bNodeRight != nullptr) {
            bNodeRight->parent = root;
        }
        if (parentLeft != nullptr) {
            parentLeft->parent = bNode;
        }
        root = bNode;
    }

    void print() const {
        vector<int>* seq = this->get_sequence();
        for (int value: *seq) {
            cout << value << " ";
        }
        cout << endl;
        delete seq;
    }

    void print_solution() const {
        vector<int>* seq = this->get_sequence();
        for (int value: *seq) {
            cout << value + 1 << " ";
        }
        cout << endl;
        delete seq;
    }
};

#endif //TSPSOLVER_TOUR_TREE_HPP
