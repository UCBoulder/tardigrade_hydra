/**
******************************************************************************
* \file tardigrade_MatrixMap.h
******************************************************************************
* A C++ library for constructing maps from row-major vectors to matrices
******************************************************************************
*/

#ifndef TARDIGRADE_MATRIXMAP_H
#define TARDIGRADE_MATRIXMAP_H

#include "Eigen/Dense"

namespace tardigradeHydra {

    template <typename T, int R, int C>
    Eigen::Map<Eigen::Matrix<T, R, C, Eigen::RowMajor> > getFixedSizeMatrixMap(T *p) {
        /*!
         * Get a matrix of type T with a fixed size of RxC from the data vector
         *
         * \param *p: A pointer to the first value of the array
         */
        return Eigen::Map<Eigen::Matrix<T, R, C, Eigen::RowMajor> >(p, R, C);
    }

    template <typename T, int R>
    Eigen::Map<Eigen::Matrix<T, R, -1, Eigen::RowMajor> > getDynamicColumnSizeMatrixMap(T *p, int C) {
        /*!
         * Get a matrix of type T with a dynamic size of RxC from the data vector
         *
         * \param *p: A pointer to the first value of the array
         * \param C: The number of columns in the map
         */
        return Eigen::Map<Eigen::Matrix<T, R, -1, Eigen::RowMajor> >(p, R, C);
    }

    template <typename T>
    Eigen::Map<Eigen::Matrix<T, -1, -1, Eigen::RowMajor> > getDynamicSizeMatrixMap(T *p, int R, int C) {
        /*!
         * Get a matrix of type T with a dynamic size of RxC from the data vector
         *
         * \param *p: A pointer to the first value of the array
         * \param R: The number of rows in the map
         * \param C: The number of columns in the map
         */
        return Eigen::Map<Eigen::Matrix<T, -1, -1, Eigen::RowMajor> >(p, R, C);
    }

    template <typename T, int R, int C>
    Eigen::Map<const Eigen::Matrix<T, R, C, Eigen::RowMajor> > getFixedSizeMatrixMap(const T *p) {
        /*!
         * Get a matrix of type T with a fixed size of RxC from the data vector
         *
         * \param *p: A pointer to the first value of the array
         */
        return Eigen::Map<const Eigen::Matrix<T, R, C, Eigen::RowMajor> >(p, R, C);
    }

    template <typename T, int R>
    Eigen::Map<const Eigen::Matrix<T, R, -1, Eigen::RowMajor> > getDynamicColumnSizeMatrixMap(const T *p, int C) {
        /*!
         * Get a matrix of type T with a dynamic size of RxC from the data vector
         *
         * \param *p: A pointer to the first value of the array
         * \param C: The number of columns in the map
         */
        return Eigen::Map<const Eigen::Matrix<T, R, -1, Eigen::RowMajor> >(p, R, C);
    }

    template <typename T>
    Eigen::Map<const Eigen::Matrix<T, -1, -1, Eigen::RowMajor> > getDynamicSizeMatrixMap(const T *p, int R, int C) {
        /*!
         * Get a matrix of type T with a dynamic size of RxC from the data vector
         *
         * \param *p: A pointer to the first value of the array
         * \param R: The number of rows in the map
         * \param C: The number of columns in the map
         */
        return Eigen::Map<const Eigen::Matrix<T, -1, -1, Eigen::RowMajor> >(p, R, C);
    }

    template <typename T, int R>
    Eigen::Map<Eigen::Vector<T, R> > getFixedSizeVectorMap(T *p) {
        /*!
         * Get a vector of type T with a fixed size of R from the data vector
         *
         * \param *p: A pointer to the first value of the array
         */

        return Eigen::Map<Eigen::Vector<T, R> >(p, R);
    }

    template <typename T, int R>
    Eigen::Map<const Eigen::Vector<T, R> > getFixedSizeVectorMap(const T *p) {
        /*!
         * Get a vector of type T with a fixed size of R from the data vector
         *
         * \param *p: A pointer to the first value of the array
         */

        return Eigen::Map<const Eigen::Vector<T, R> >(p, R);
    }

    template <typename T>
    Eigen::Map<Eigen::Vector<T, -1> > getDynamicSizeVectorMap(T *p, int R) {
        /*!
         * Get a vector of type T with a fixed size of R from the data vector
         *
         * \param *p: A pointer to the first value of the array
         * \param R: The number of rows in the map
         */

        return Eigen::Map<Eigen::Vector<T, -1> >(p, R);
    }

    template <typename T>
    Eigen::Map<const Eigen::Vector<T, -1> > getDynamicSizeVectorMap(const T *p, int R) {
        /*!
         * Get a vector of type T with a fixed size of R from the data vector
         *
         * \param *p: A pointer to the first value of the array
         * \param R: The number of rows in the map
         */

        return Eigen::Map<const Eigen::Vector<T, -1> >(p, R);
    }

}  // namespace tardigradeHydra

#include "tardigrade_MatrixMap.cpp"

#endif
