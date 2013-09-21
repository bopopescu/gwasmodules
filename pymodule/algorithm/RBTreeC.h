/*
 * 2013.09.18 Yu Huang, copyright. a red-black tree template, c++
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

#define COLOR(color) (color == 0) ? "red" : "black"

namespace YHuang {

static const short BLACK_ = 1;
static const short RED_ = 0;
static const short LEFT_ = 100;
static const short RIGHT_ = 200;

template <typename keyType, typename dataType>
class RBTreeNode {
private:
	RBTreeNode<keyType, dataType> *parent_;
	RBTreeNode<keyType, dataType> *left_;
	RBTreeNode<keyType, dataType> *right_;
	keyType* keyPtr;
	dataType* dataPtr;
	unsigned short color_;
public:
	RBTreeNode(RBTreeNode<keyType, dataType>* _parent, keyType* _keyPtr, dataType* _dataPtr): parent_(_parent),
		keyPtr(_keyPtr), dataPtr(_dataPtr){
		/*
		 * key_, data_ are references, have to be initialized in the away above.
		 */
		this->left_ = NULL;
		this->right_ = NULL;
		this->color_ = RED_;
	}
	~RBTreeNode(){
		delete parent_;
		delete this->left_;
		delete this->right_;
	}

	RBTreeNode<keyType, dataType> *getParent(){
		return this->parent_;
	}
	RBTreeNode<keyType, dataType> *getLeft(){
		return this->left_;
	}
	RBTreeNode<keyType, dataType> *getRight(){
		return this->right_;
	}
	dataType* getData(){
		return this->dataPtr;
	}
	keyType* getKey(){
		return this->keyPtr;
	}
	short getColor(){
		return this->color_;
	}
	void setLeft(RBTreeNode<keyType, dataType> *nodePtr){
		this->left_ = nodePtr;
	}
	void setRight(RBTreeNode<keyType, dataType> *nodePtr){
		this->right_ = nodePtr;
	}
	void setParent(RBTreeNode<keyType, dataType> *nodePtr){
		this->parent_ = nodePtr;
	}
	void setData(dataType* dataPtr){
		this->dataPtr = dataPtr;
	}
	void setKey(keyType* keyPtr){
		this->keyPtr = keyPtr;
	}

	void setColor(short color){
		this->color_ = color;

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
			//printf("Rotate Right Of %d\n", nodePtr->keyPtr);
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
			//printf("Rotate Left Of %d\n", nodePtr->keyPtr);
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
		RBTreeNode<keyType, dataType> *treeInsert_(RBTreeNode<keyType, dataType> *nodePtr, keyType* keyPtr, dataType* dataPtr){
			RBTreeNode<keyType, dataType> *inserted;
			if (nodePtr==NULL){	//2013.09.20 YH. first time insert, root_ is NULL
				//inserted = createNode_(nodePtr, LEFT_, keyPtr, dataPtr);
				root_ = new RBTreeNode<keyType, dataType>(NULL, keyPtr, dataPtr);
				inserted = root_;
				insertFix_(inserted);
			}
			else if (*keyPtr < *nodePtr->getKey()) {
				if (nodePtr->getLeft() == NULL) {
					inserted = createNode_(nodePtr, LEFT_, keyPtr, dataPtr);
					insertFix_(inserted);
				} else {
					inserted = treeInsert_(nodePtr->getLeft(), keyPtr, dataPtr);
				}
			} else {
				if (nodePtr->getRight() == NULL) {
					inserted = createNode_(nodePtr, RIGHT_, keyPtr, dataPtr);
					//printf("Inserted %d On Right\n", keyPtr);
					insertFix_(inserted);
				} else {
					inserted = treeInsert_(nodePtr->getRight(), keyPtr, dataPtr);
				}
			}
			//FIX Location
			return inserted;
		}
		RBTreeNode<keyType, dataType> *createNode_(RBTreeNode<keyType, dataType> *parent, short loc, keyType* keyPtr, dataType* dataPtr){
			RBTreeNode<keyType, dataType> *tmp = new RBTreeNode<keyType, dataType>(parent, keyPtr, dataPtr);
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
							continue;
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
				} else {
					break;
				}
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
		RBTreeNode<keyType, dataType> *queryTreeRecur_(RBTreeNode<keyType, dataType> *nodePtr, keyType* keyPtr){
			if (nodePtr == NULL) {
				return NULL;
			}
			if (*keyPtr < *nodePtr->getKey()) {
				return queryTreeRecur_(nodePtr->getLeft(), keyPtr);
			} else if (*keyPtr > *nodePtr->getKey()) {
				return queryTreeRecur_(nodePtr->getRight(), keyPtr);
			} else {
				return nodePtr;
			}
		}
		void printTreeRecur_(RBTreeNode<keyType, dataType> *nodePtr){
			if (nodePtr == NULL) {
				return;
			}
			printTreeRecur_(nodePtr->getLeft());
			std::cout<<boost::format("%1% ")%  nodePtr->getKey();
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
	public:
		//Basic Functions
		RBTree(){
			root_ = NULL;

		}
		RBTree(keyType* keyPtr, dataType* dataPtr){
			root_ = new RBTreeNode<keyType, dataType>(NULL, keyPtr, dataPtr);
		}
		~RBTree(){
			delete this->root_;
		}
		//rbtreePtr createRBTree(int, dataType );
		//void destroyRBTree();
		RBTreeNode<keyType, dataType> *insertNode(keyType* keyPtr, dataType* dataPtr){
			return treeInsert_(this->root_, keyPtr, dataPtr);

		}
		long deleteNodeByKey(keyType* keyPtr){
			RBTreeNode<keyType, dataType> *nodePtr = queryTree(keyPtr);
			return deleteNode(nodePtr);
		}
		long deleteNode(RBTreeNode<keyType, dataType> *nodePtr){
			//RBTreeNode<keyType, dataType> *nodePtr = queryTree(keyPtr);
			RBTreeNode<keyType, dataType> *y;
			RBTreeNode<keyType, dataType> *x;
			if (nodePtr == NULL)
				return 0;
			if (nodePtr->getLeft() == NULL || nodePtr->getRight() == NULL) {
				y = nodePtr;
			} else {
				y = getSuccessor_(nodePtr);
			}
			if (y->getLeft() != NULL) {
				x = y->getLeft();
			} else {
				x = y->getRight();
			}
			if (x != NULL) {	//remove y from chain
				x->setParent(y->getParent());
			}
			if (y->getParent() == NULL) {
				this->root_ = x;
			} else {
				if (y == y->getParent()->getLeft()) {
					y->getParent()->setLeft(x);
				} else {
					y->getParent()->setRight(x);
				}
			}
			if (y != nodePtr) {
				nodePtr->setKey(y->getKey());
				nodePtr->setData(y->getData());
			}
			delete y;
			return 0;
		}
		short isValidRedBlackTree(){
			return isValidRedBlackTreeRecur_(this->root_);
		}
		RBTreeNode<keyType, dataType> *queryTree(keyType* keyPtr){
			return queryTreeRecur_(this->root_, keyPtr);
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
			printTreeRecur_(this->root_);
		}
		void printPaths(){
			RBTreeNode<keyType, dataType> *path[1000];
			printPathsRecur_(this->root_, path, 0);
		}


};

}

#endif /* RBTREES_H_ */
