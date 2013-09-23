/*
 * 2013.09.18 Yu Huang, copyright. a red-black tree template, c++
 * still contains bugs, use RedBlackTree.h instead.
 *
 * Built on top of asoliman's work which was created on: Feb 1, 2009 https://code.google.com/p/rbtrees/.
 * 	his version contains bugs. i.e. in deleteNode() and sibling_() is not implemented.
 */

#ifndef RBTREES_H_
#define RBTREES_H_
#include <inttypes.h>
#include <stdio.h>
#include <iostream>
#include <boost/format.hpp>
#include <ostream>
using namespace std;

#define COLOR(color) (color == 0) ? "red" : "black"

namespace YHuang {

static const short RED_ = 0;
static const short BLACK_ = 1;
static const short LEFT_ = 100;
static const short RIGHT_ = 200;

template <typename keyType, typename dataType>
class RBTreeNode {
	// for output
	friend ostream& operator<<(ostream& out, RBTreeNode& node){
			//out << boost::format(" RBTreeNode object key=%1%, dataPtr=%2%, left=%3%, right=%4%, parent=%5%, color=%6% ")%
			//		node.getKey() % node.getDataPtr() % node.getLeft() %
			//		node.getRight() % node.getParent() % node.getColor();
			out << boost::format(" RBTreeNode object");
			return out;
		}
public:
	RBTreeNode<keyType, dataType> *parent;
	RBTreeNode<keyType, dataType> *left;
	RBTreeNode<keyType, dataType> *right;
	keyType key;
	dataType* dataPtr;
	unsigned short color;
	RBTreeNode(RBTreeNode<keyType, dataType>* _parent, keyType _key, dataType* _dataPtr): parent(_parent),
		key(_key), dataPtr(_dataPtr){
		/*
		 */
		this->left = NULL;
		this->right = NULL;
		this->color = RED_;
	}
	~RBTreeNode(){
		free(parent);	//if it's NULL, nothing will happen
		free(this->left);
		free(this->right);
	}

	RBTreeNode<keyType, dataType> *getParent(){
		return this->parent;
	}
	RBTreeNode<keyType, dataType> *getLeft(){
		return this->left;
	}
	RBTreeNode<keyType, dataType> *getRight(){
		return this->right;
	}
	dataType* getDataPtr(){
		return this->dataPtr;
	}
	keyType getKey(){
		return this->key;
	}
	short getColor(){
		return this->color;
	}
	void setLeft(RBTreeNode<keyType, dataType> *nodePtr){
		this->left = nodePtr;
	}
	void setRight(RBTreeNode<keyType, dataType> *nodePtr){
		this->right = nodePtr;
	}
	void setParent(RBTreeNode<keyType, dataType> *nodePtr){
		this->parent = nodePtr;
	}
	void setDataPtr(dataType* dataPtr){
		this->dataPtr = dataPtr;
	}
	void setKey(keyType key){
		this->key = key;
	}

	void setColor(short color){
		this->color = color;

	}


};

template <typename keyType, typename dataType>
class RBTree {
	private:
		// Helper Functions
		RBTreeNode<keyType, dataType> *root_;
		RBTreeNode<keyType, dataType> *getMinimum_(RBTreeNode<keyType, dataType> *nodePtr){
			if (nodePtr == NULL) {
				return NULL;
			}
			if (nodePtr->getLeft() != NULL) {
				return getMinimum_(nodePtr->getLeft());
			}
			return nodePtr;
		}
		RBTreeNode<keyType, dataType> *getMaximum_(RBTreeNode<keyType, dataType> *nodePtr){
			if (nodePtr == NULL) {
				return NULL;
			}
			if (nodePtr->getRight() != NULL) {
				return getMaximum_(nodePtr->getRight());
			}
			return nodePtr;
		}
		RBTreeNode<keyType, dataType> *grandparent_(RBTreeNode<keyType, dataType> *nodePtr){
			if ((nodePtr != NULL) && (nodePtr->getParent() != NULL)) {
				return nodePtr->getParent()->getParent();
			} else {
				return NULL;
			}
		}
		RBTreeNode<keyType, dataType> *uncle_(RBTreeNode<keyType, dataType> *nodePtr){
			RBTreeNode<keyType, dataType> *myGPNode = grandparent_(nodePtr);
			if (myGPNode == NULL) {
				return NULL;
			}
			if (nodePtr->getParent() == myGPNode->getLeft()) {
				return myGPNode->getRight();
			} else {
				return myGPNode->getLeft();
			}
		}
		RBTreeNode<keyType, dataType> *getSuccessor_(RBTreeNode<keyType, dataType> *nodePtr){
			/*
			 * 2013.09.21 successor is the node that is the next bigger one after nodePtr
			 */
			if (nodePtr->getRight() != NULL) {
				//leftmost of the right child
				return getMinimum_(nodePtr->getRight());
			}
			//first right-side parent of all left-side ancestors
			RBTreeNode<keyType, dataType> *x = nodePtr;
			RBTreeNode<keyType, dataType> *y = nodePtr->getParent();
			while (y != NULL && x == y->getRight()) {
				x = y;
				y = y->getParent();
			}
			return y;
		}
		RBTreeNode<keyType, dataType> *getPredecessor_(RBTreeNode<keyType, dataType> *nodePtr){
			if (nodePtr->getLeft() != NULL) {
				return getMaximum_(nodePtr->getLeft());
			}
			RBTreeNode<keyType, dataType> *x = nodePtr;
			RBTreeNode<keyType, dataType> *y = nodePtr->getParent();
			while (y != NULL && x == y->getLeft()) {
				x = y;
				y = y->getParent();
			}
			return y;
		}
		RBTreeNode<keyType, dataType> *sibling_(RBTreeNode<keyType, dataType> *nodePtr){
			/*
			 * 2013.09.21 YH
			 */
			if (nodePtr->getParent()!=NULL){
				if (nodePtr == nodePtr->getParent()->getLeft())
					return nodePtr->getParent()->getRight();
				else
					return nodePtr->getParent()->getLeft();
			}
			else{
				return NULL;
			}
		}

		short isValidRedBlackTreeRecur_(RBTreeNode<keyType, dataType> *nodePtr){
			if (nodePtr == NULL)
				return 1;
			if (!isValidRedBlackTreeRecur_(nodePtr->getLeft())) {
				return 0;
			}
			if (!isValidRedBlackTreeRecur_(nodePtr->getRight())) {
				return 0;
			}
			if (nodePtr->getParent() == NULL && nodePtr->getColor() == BLACK_) {
				return 1;
			}
			if (nodePtr->getColor() == RED_ && nodePtr->getParent() != NULL && nodePtr->getParent()->getColor()
					== BLACK_) {
				return 1;
			}
			if (nodePtr->getColor() == BLACK_ && nodePtr->getParent() != NULL && nodePtr->getParent()->getColor()
					== BLACK_) {
				return 1;
			}
			return 0;
		}
		short rotateRight_(RBTreeNode<keyType, dataType> *nodePtr){
			//printf("Rotate Right Of %d\n", nodePtr->key);
			//cool operation...
			if (nodePtr->getLeft() == NULL)
				return 1;
			RBTreeNode<keyType, dataType> *leftNode = nodePtr->getLeft();
			RBTreeNode<keyType, dataType> *correctParent = nodePtr->getParent();
			//let's fix the parent's links first...
			if (correctParent != NULL) {
				if (nodePtr == correctParent->getLeft()) {
					correctParent->setLeft(leftNode);
				} else {
					correctParent->setRight(leftNode);
				}
			}
			nodePtr->setLeft(leftNode->getRight());
			leftNode->setRight(nodePtr);
			leftNode->setParent(nodePtr->getParent());
			nodePtr->setParent(leftNode);

			if (this->root_ == nodePtr) {
				this->root_ = leftNode;
				leftNode->setParent(NULL);
			}
			return 0;
		}
		short rotateLeft_(RBTreeNode<keyType, dataType> *nodePtr){
			//printf("Rotate Left Of %d\n", nodePtr->key);
			//cool operation...
			if (nodePtr->getRight() == NULL)
				return 1;
			RBTreeNode<keyType, dataType> *rightNode = nodePtr->getRight();
			RBTreeNode<keyType, dataType> *correctParent = nodePtr->getParent();
			//let's fix the parent's links first...
			if (correctParent != NULL) {
				if (nodePtr == correctParent->getLeft()) {
					correctParent->setLeft(rightNode);
				} else {
					correctParent->setRight(rightNode);
				}
			}
			nodePtr->setRight(rightNode->getLeft());
			rightNode->setLeft(nodePtr);
			rightNode->setParent(nodePtr->getParent());
			nodePtr->setParent(rightNode);
			if (this->root_ == nodePtr) {
				this->root_ = rightNode;
				rightNode->setParent(NULL);
			}
			return 0;
		}
		short isLeaf_(RBTreeNode<keyType, dataType> *nodePtr){
			if (nodePtr->getLeft() == NULL && nodePtr->getRight() == NULL)
				return 1;
			return 0;

		}
		RBTreeNode<keyType, dataType> *treeInsert_(RBTreeNode<keyType, dataType> *parentNodePtr, keyType key, dataType* dataPtr){
			RBTreeNode<keyType, dataType> *inserted=NULL;
			if (parentNodePtr==NULL){	//2013.09.20 YH. first time insert, root_ is NULL
				inserted = new RBTreeNode<keyType, dataType>(NULL, key, dataPtr);
				//inserted = createNode_(parentNodePtr, LEFT_, key, dataPtr);
				//this->root_ = new RBTreeNode<keyType, dataType>(NULL, key, dataPtr);
				this->root_ = inserted;
				//this->root_->setColor(BLACK_);
				//insertFix_(inserted);

			}
			else if (key < parentNodePtr->getKey()) {
				if (parentNodePtr->getLeft() == NULL) {
					inserted = createNode_(parentNodePtr, LEFT_, key, dataPtr);
					insertFix_(inserted);
				} else {
					inserted = treeInsert_(parentNodePtr->getLeft(), key, dataPtr);
				}
			} else {
				if (parentNodePtr->getRight() == NULL) {
					inserted = createNode_(parentNodePtr, RIGHT_, key, dataPtr);
					//printf("Inserted %d On Right\n", key);
					insertFix_(inserted);
				} else {
					inserted = treeInsert_(parentNodePtr->getRight(), key, dataPtr);
				}
			}
			//FIX Location
			return inserted;
		}
		RBTreeNode<keyType, dataType> *createNode_(RBTreeNode<keyType, dataType> *parent, short loc, keyType key, dataType* dataPtr){
			RBTreeNode<keyType, dataType> *tmp = new RBTreeNode<keyType, dataType>(parent, key, dataPtr);
			if (loc == LEFT_)
				parent->setLeft(tmp);
			else
				parent->setRight(tmp);
			return tmp;
		}
		void insertFix_(RBTreeNode<keyType, dataType> *nodePtr){
			RBTreeNode<keyType, dataType> *z = nodePtr;
			RBTreeNode<keyType, dataType> *y = NULL;
			while (z != NULL && z->getParent() != NULL && z->getParent()->getColor() == RED_) {
				if (z->getParent()->getParent() != NULL) {
					if (z->getParent() == z->getParent()->getParent()->getLeft()) {
						//case 1,2,3
						y = uncle_(z);
						if (y != NULL && y->getColor() == RED_) {
							//case 1
							//RECOLOR
							z->getParent()->setColor(BLACK_);
							y->setColor(BLACK_);
							z->getParent()->getParent()->setColor(RED_);
							z = z->getParent()->getParent();
							//continue;
						} else {
							if (z == z->getParent()->getRight()) {
								//case 2
								z = z->getParent();
								rotateLeft_(z);
							}
							z->getParent()->setColor(BLACK_);
							z->getParent()->getParent()->setColor(RED_);
							rotateRight_(z->getParent()->getParent());
							//case 3
						}
					} else {
						//case 4,5,6
						y = uncle_(z);

						if (y != NULL && y->getColor() == RED_) {
							//case 4
							//RECOLOR
							z->getParent()->setColor(BLACK_);
							y->setColor(BLACK_);
							z->getParent()->getParent()->setColor(RED_);
							z = z->getParent()->getParent();
						} else {
							if (z == z->getParent()->getLeft()) {
								//case 5
								z = z->getParent();
								rotateRight_(z);
							}
							//case 6
							z->getParent()->setColor(BLACK_);
							z->getParent()->getParent()->setColor(RED_);
							rotateLeft_(z->getParent()->getParent());
						}
					}
				}
				//else {
				//	break;
				//}
			}
			this->root_->setColor(BLACK_);
		}
		long count_(RBTreeNode<keyType, dataType> *nodePtr, long num){
			if (nodePtr == NULL) {
				return num;
			}
			return count_(nodePtr->getLeft(), count_(nodePtr->getRight(), ++num));
		}
		long maxDepthRecur_(RBTreeNode<keyType, dataType> *nodePtr){
			if (nodePtr == NULL) {
				return 0;
			}
			long leftDepth = maxDepthRecur_(nodePtr->getLeft());
			long rightDepth = maxDepthRecur_(nodePtr->getRight());
			if (leftDepth > rightDepth) {
				return leftDepth + 1;
			} else {
				return rightDepth + 1;
			}
		}
		RBTreeNode<keyType, dataType> *queryTreeRecur_(RBTreeNode<keyType, dataType> *nodePtr, keyType key){
			if (nodePtr == NULL) {
				return NULL;
			}
			if (key < nodePtr->key) {
				return queryTreeRecur_(nodePtr->getLeft(), key);
			} else if (key > nodePtr->key) {
				return queryTreeRecur_(nodePtr->getRight(), key);
			} else {
				return nodePtr;
			}
		}
		void printTreeRecur_(RBTreeNode<keyType, dataType> *nodePtr){
			if (nodePtr == NULL) {
				return;
			}
			printTreeRecur_(nodePtr->getLeft());
			std::cout<<boost::format("%1% ") % nodePtr->key;
			printTreeRecur_(nodePtr->getRight());
		}
		void printPathsRecur_(RBTreeNode<keyType, dataType> *nodePtr, RBTreeNode<keyType, dataType>** path, int pathLen){
			if (nodePtr == NULL)
					return;

			path[pathLen++] = nodePtr;
			if (isLeaf_(nodePtr)) {
				long i;
				for (i = 0; i < pathLen; i++) {
					std::cout<<boost::format("%1% ")% path[i]->getKey();
				}
				printf("\n");
				return;
			}
			printPathsRecur_(nodePtr->getLeft(), path, pathLen);
			printPathsRecur_(nodePtr->getRight(), path, pathLen);

		}
		void deleteFixup_(RBTreeNode<keyType, dataType> *x) {
			/*************************************
			 *  maintain Red-Black tree balance  *
			 *  after
			 *  1. node x's parent y, has been moved into the deleted node's place
			 *  2. x has replaced y's place
			 *  3. y's color is black
			 *  4. x's color is also black
			 *  5. x is not root
			 *************************************/
			while (x!=NULL && x != root_ && x->getColor() == BLACK_) {
				if (x == x->getParent()->getLeft()) {
					RBTreeNode<keyType, dataType> *siblingPtr = x->getParent()->getRight();	//get sibling
					//sibling must have valid (non-leaf) left and right because both x's parent and x is black, (property 5)
					if (siblingPtr->getColor() == RED_) {	//case 2
						siblingPtr->setColor(BLACK_);
						x->getParent()->setColor(RED_);
						rotateLeft_(x->getParent());
						siblingPtr = x->getParent()->getRight();
					}
					if (siblingPtr->getLeft()->getColor() == BLACK_ && siblingPtr->getRight()->getColor() == BLACK_) {
						//meet property 5 to maintain the same number of black nodes from x's parent
						siblingPtr->setColor(RED_);
						x = x->getParent();
					} else {
						if (siblingPtr->getRight()->getColor() == BLACK_) {
							siblingPtr->getLeft()->setColor(BLACK_);
							siblingPtr->setColor(RED_);
							rotateRight_(siblingPtr);
							siblingPtr = x->getParent()->getRight();
						}
						siblingPtr->setColor(x->getParent()->getColor());
						x->getParent()->setColor(BLACK_);
						siblingPtr->getRight()->setColor(BLACK_);
						rotateLeft_(x->getParent());
						x = root_;
					}
				} else {
					RBTreeNode<keyType, dataType> *siblingPtr = x->getParent()->getLeft();	//get sibling
					if (siblingPtr==NULL){
						std::cerr << boost::format("deleteFixup_() Error: siblingPtr (address=%1%, %2%) of x (addr=%3%, %4%) >> is NULL, but can't happen as x is black.")%
								siblingPtr % (*siblingPtr) % x % (*x) << std::endl;
					}
					if (siblingPtr->getColor() == RED_) {
						siblingPtr->setColor(BLACK_);
						x->getParent()->setColor(RED_);
						rotateRight_(x->getParent());
						siblingPtr = x->getParent()->getLeft();
					}
					if (siblingPtr->getRight()->getColor() == BLACK_ && siblingPtr->getLeft()->getColor() == BLACK_) {
						siblingPtr->setColor(RED_);
						x = x->getParent();
					} else {
						if (siblingPtr->getLeft()->getColor() == BLACK_) {
							siblingPtr->getRight()->setColor(BLACK_);
							siblingPtr->setColor(RED_);
							rotateLeft_(siblingPtr);
							siblingPtr = x->getParent()->getLeft();
						}
						siblingPtr->setColor(x->getParent()->getColor());
						x->getParent()->setColor(BLACK_);
						siblingPtr->getLeft()->setColor(BLACK_);
						rotateRight_(x->getParent());
						x = root_;
					}
				}
			}
			if (x!=NULL){
				x->setColor(BLACK_);
			}
		}

	public:
		//Basic Functions
		RBTree(){
			root_ = NULL;

		}
		RBTree(keyType key, dataType* dataPtr){
			root_ = new RBTreeNode<keyType, dataType>(NULL, key, dataPtr);
		}
		~RBTree(){
			delete this->root_;
		}
		//rbtreePtr createRBTree(int, dataType* );
		//void destroyRBTree();
		RBTreeNode<keyType, dataType> *insertNode(keyType key, dataType* dataPtr){
			RBTreeNode<keyType, dataType>* inserted = new RBTreeNode<keyType, dataType>(NULL, key, dataPtr);
			RBTreeNode<keyType, dataType>* y = NULL;
			RBTreeNode<keyType, dataType>* x = root_;
			while (x!=NULL){
				y=x;
				if (inserted->key < x->key){
					x = x->getLeft();
				}
				else{
					x = x->getRight();
				}

			}
			inserted->setParent(y);
			if (y==NULL){
				root_ = inserted;
			}else if (inserted->getKey() < y->getKey()){
				y->setLeft(inserted);

			}else{
				y->setRight(inserted);
			}
			inserted->setLeft(NULL);
			inserted->setRight(NULL);
			inserted->setColor(RED_);
			insertFix_(inserted);
			//return treeInsert_(this->root_, key, dataPtr);
			return inserted;

		}
		long deleteNodeByKey(keyType key){
			RBTreeNode<keyType, dataType> *nodePtr = queryTree(key);
			return deleteNode(nodePtr);
		}

		long deleteNode(RBTreeNode<keyType, dataType> *nodePtr){
			//RBTreeNode<keyType, dataType> *nodePtr = queryTree(key);
			RBTreeNode<keyType, dataType> *y;
			RBTreeNode<keyType, dataType> *x;
			if (nodePtr == NULL)
				return 0;
			if (nodePtr->getLeft() == NULL || nodePtr->getRight() == NULL) {
				y = nodePtr;
			} else {
				y = getSuccessor_(nodePtr);
			}
			if (y->getLeft() != NULL) {	//try to get a non-leaf child, but could still be a leaf (NULL) if both are NULL
				x = y->getLeft();
			} else {
				x = y->getRight();
			}
			if (x != NULL) {	//remove y from chain, replaces it with x
				x->setParent(y->getParent());
			}
			if (y->getParent() == NULL) {
				this->root_ = x;
			} else {
				if (y == y->getParent()->getLeft()) {
					y->getParent()->setLeft(x);
				} else if (y==y->getParent()->getRight()) {
					y->getParent()->setRight(x);
				} else{
					std::cerr << boost::format("deleteNode() Error: node y, << address=%1%, %2%>> , is neither left , nor right of its parent << address=%3%, %4% >>.\n") % y % (*y) % y->getParent() % (*y->getParent()) << std::endl;
				}
			}
			if (y != nodePtr) {	//put y's key and dataPtr into nodePtr, since y is next in line to succeed nodePtr, no change.
				nodePtr->setKey(y->getKey());
				nodePtr->setDataPtr(y->getDataPtr());
			}
			if (y->getColor() == BLACK_){	//x (y's child) has replaced y and replacing a black node (y) requires rebalancing
				deleteFixup_(x);
			}
			//gets rid of y
			y->setParent(NULL);
			y->setLeft(NULL);
			y->setRight(NULL);
			y=NULL;
			free(y);
			return 0;
		}
		short isValidRedBlackTree(){
			return isValidRedBlackTreeRecur_(this->root_);
		}
		RBTreeNode<keyType, dataType> *queryTree(keyType key){
			return queryTreeRecur_(this->root_, key);
		}
		long size(){
			return count_(this->root_, 0);
		}
		long maxDepth(){
			return maxDepthRecur_(this->root_);
		}
		RBTreeNode<keyType, dataType> *getMinimum(){
			return getMinimum_(this->root_);

		}
		RBTreeNode<keyType, dataType> *getMaximum(){
			return getMaximum_(this->root_);

		}

		//Visualization Functions
		void printTree(){
			if (root_!=NULL)
				printTreeRecur_(this->root_);
		}
		void printPaths(){
			RBTreeNode<keyType, dataType> *path[1000];
			printPathsRecur_(this->root_, path, 0);
		}

		bool isNULLNode(RBTreeNode<keyType, dataType> * nodePtr){
			return nodePtr==NULL;
		}
};

}

#endif /* RBTREES_H_ */
