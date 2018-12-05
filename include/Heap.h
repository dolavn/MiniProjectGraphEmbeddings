#include <vector>
#include <string>
#include <iostream>

#ifndef HEAP_H_
#define HEAP_H_
/**
An heap data structure.

@author Dolav Nitay
*/
template <class T>
class Heap
{
public:
	/**
	Creates a new empty heap object.

	@param cmp A function to compare between two T objects.
	@param print A function to print each T object.
	*/
	Heap(bool(*cmp)(T, T), std::string(*print)(T));
	/**
	Receives a vector of T objects, and creates a heap data structure form it in O(n) time.

	@param arr The vector from which a heap would be created
	@param cmp A function to compare between two T objects
	@param print A function to print each T object
	*/
	Heap(std::vector<T> arr, bool(*cmp)(T, T), std::string(*print)(T));

	/**
	Receives an array of T objects, which already satisfies heap conditions. 
	Forms a heap from it in O(1) time.

	@param arr The vector from which a heap would be created
	@param cmp A function to compare between two T objects
	@param print A function to print each T object
	@param sorted Indicates whether the given vector satisfies the heap constraints
	*/
	Heap(std::vector<T> heap, bool(*cmp)(T, T), std::string(*print)(T),bool sorted);
	/**
	Creates a deep copy of a given heap.

	@param other The other heap.
	*/
	Heap(const Heap<T>& other);

	/**
	Destructs this heap object, frees the memory allocated by it.
	*/
	~Heap();
	/**
	Removes the head of the heap and returns it to the user in O(log n) time.

	@return The head of the heap
	*/
	T popHead();
	/**
	Inserts a new element to the heap in O(log n) time.

	@param value The value to be inserted
	@return A pointer to the index of the newly inserted element in the heap.
	*/
	int* pushHeap(T value);
	/**
	Returns the size of the heap.

	@return The size of the heap
	*/
	unsigned int getSize() { return vec.size(); }
	/**
	Updates the value of given object in the heap in O(log n) time.

	@param ind The current index in the heap of the object.
	@param newVal The new value of the object in the heap
	*/
	void setValue(int ind, T newVal);
	/**
	Prints the heap.
	*/
	void printHeap();
	/**
	Used when constructing an heap from a given vector. Returns a vector of pointers to the current indexes in the heap
	of each object.

	@return A vector of the pointers to the indexes in the heap of each object in the original array.
	*/
	std::vector<int*>& getPointerVec() { return pointersVec; }
private:
	bool(*cmp)(T, T); //Function to compare between heap elements
	std::string(*print)(T); //Function to print heap elements
	std::vector<std::pair<int*, T>> vec; //The heap
	std::vector<int*> pointersVec; //The pointers vec

	T getElem(int ind); //Returns the elem in index ind in the heap
	int getFather(int i) { return (i - 1) / 2; } //Returns the index of the father of a given index in the heap.
	int getLeftSon(int i) { return (2 * i) + 1; } //Returns the index of the left child of a given index in the heap.
	int getRightSon(int i) { return (2 * i) + 2; } //Returns the index of the right child of a given index in the heap.
	void downHeapify(int ind); //Pushes the object in index ind downwards in the heap, until it satisfies the heap requirements
	void upHeapify(int ind); //Pushes the object in index ind upwards in the heap, until it satisfies the heap requirements
	void buildHeap(std::vector<T> arr); //Builds the heap from a given vector of objects in O(n) time
	void swap(int ind1, int ind2); //Swaps between two values in the heap
};



template<class T>
Heap<T>::Heap(bool(*cmp)(T, T), std::string(*print)(T)) :cmp(cmp), print(print), vec(), pointersVec() {
}

template<class T>
Heap<T>::Heap(std::vector<T> arr, bool(*cmp)(T, T), std::string(*print)(T)) : cmp(cmp), print(print), vec(arr.size()), pointersVec(arr.size()) {
	buildHeap(arr);
}

template<class T>
Heap<T>::Heap(std::vector<T> heap, bool(*cmp)(T, T), std::string(*print)(T), bool sorted) : cmp(cmp), print(print), vec(heap.size()), pointersVec(heap.size()) {
	for (unsigned int i = 0; i<vec.size(); i++) { //O(n)
		int* index = new int(i);
		std::pair<int*, T> currPair = std::pair<int*, T>(index, heap[i]);
		vec[i] = currPair;
		pointersVec[i] = index;
	}
}

template<class T>
Heap<T>::Heap(const Heap<T>& other) : cmp(other.cmp), print(other.print), vec(other.vec.size()), pointersVec(other.vec.size()) {
	for (unsigned int i = 0; i<vec.size(); i++) {
		std::pair<int*, T> currPair;
		currPair.first = new int(i);
		currPair.second = other.vec[i].second;
		vec[i] = currPair;
		pointersVec[i] = currPair.first;
	}
}

template<class T>
Heap<T>::~Heap() {
	for (unsigned int i = 0; i<vec.size(); i++) {
		std::pair<int*, T> currPair = vec[i];
		if (nullptr != currPair.first) {
			delete(currPair.first);
			currPair.first = nullptr;
		}
	}
}

template<class T>
T Heap<T>::popHead() {
	T ans = getElem(0);
	swap(0, vec.size() - 1);
	int*& removedIndex = vec[vec.size() - 1].first;
	if (removedIndex != nullptr) {delete(removedIndex); removedIndex = nullptr; }
	vec.pop_back();
	if (vec.size() > 0) {
		downHeapify(0);
	}
	return ans;
}

template<class T>
int* Heap<T>::pushHeap(T value) {
	int* ans = new int(vec.size());
	std::pair<int*, T> newPair(ans,value);
	vec.push_back(newPair);
	upHeapify(vec.size() - 1);
	return ans;
}

template<class T>
void Heap<T>::setValue(int ind, T newVal) {
	vec[ind].second = newVal;
	upHeapify(ind);
}

template<class T>
void Heap<T>::downHeapify(int ind) {
	int minSon = ind; //The index of the son with the minimum value.
	unsigned int sonL = getLeftSon(ind);
	unsigned int sonR = getRightSon(ind);
	if (sonL<vec.size() && cmp(getElem(sonL), getElem(minSon))) { minSon = sonL; }
	if (sonR<vec.size() && cmp(getElem(sonR), getElem(minSon))) { minSon = sonR; }
	while (cmp(getElem(minSon), getElem(ind))) { //ind should be before min son
		swap(minSon, ind); //swaps between ind and the minimal son
		ind = minSon; //Updates indexes
		minSon = ind;
		sonL = getLeftSon(ind);
		sonR = getRightSon(ind);
		if (sonL<vec.size() && cmp(getElem(sonL), getElem(minSon))) { minSon = sonL; }
		if (sonR<vec.size() && cmp(getElem(sonR), getElem(minSon))) { minSon = sonR; }
	}

}

template<class T>
void Heap<T>::upHeapify(int ind) {
	int father = getFather(ind);
	while (ind>0 && cmp(getElem(ind), getElem(father))) {
		swap(ind, father);
		ind = father;
		father = getFather(ind);
		if (father<0) { throw 1; }
	}
}

template<class T>
void Heap<T>::buildHeap(std::vector<T> arr) {
	for (unsigned int i = 0; i<vec.size(); i++) { //O(n)
		int* index = new int(i);
		std::pair<int*, T> currPair = std::pair<int*, T>(index, arr[i]);
		vec[i] = currPair;
		pointersVec[i] = index;
	}
	for (int i = (arr.size() / 2); i >= 0; i--) {
		downHeapify(i);
	}
}

template<class T>
void Heap<T>::swap(int ind1, int ind2) {
	std::pair<int*, T> tmp = vec[ind1];
	vec[ind1] = vec[ind2];
	vec[ind2] = tmp;
	*(vec[ind1].first) = ind1; //Updates the indexes in the pointers
	*(vec[ind2].first) = ind2;
}

template<class T>
T Heap<T>::getElem(int ind) {
	return vec[ind].second;
}

template<class T>
void Heap<T>::printHeap() {
	for (unsigned int i = 0; i < pointersVec.size(); i++) {
		std::cout << *pointersVec[i] << std::endl;
	}
	for (unsigned int i = 0; i<vec.size(); i++) {
		std::cout << "i:"; std::cout << i;
		std::cout << " ind:"; std::cout << *(vec[i].first);
		std::cout << " ";
		std::cout << print(getElem(i)); std::cout << std::endl;
	}
	std::cout << std::endl;
}

#endif