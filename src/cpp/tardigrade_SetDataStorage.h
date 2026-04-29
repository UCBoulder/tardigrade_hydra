/**
 ******************************************************************************
 * \file tardigrade_SetDataStorage.h
 ******************************************************************************
 * A C++ library for constructing caching data storage objects
 ******************************************************************************
 */

#ifndef TARDIGRADE_SETDATASTORAGE
#define TARDIGRADE_SETDATASTORAGE

#include "tardigrade_CoreDefinitions.h"
#include "tardigrade_MatrixMap.h"
#include "tardigrade_error_tools.h"

/*!
 * \brief Declares a named SetDataStorage getter function
 * \param getname: The name of the getter function
 * \param varname: The name of the variable
 * \param classtype: The type of the class to be defined
 * \param vartype: The ctype of the variable
 */
#define TARDIGRADE_HYDRA_DECLARE_NAMED_SETDATASTORAGE_GETTER(getname, varname, classtype, vartype) \
    classtype<vartype> getname() {                                                                 \
        /*!                                                                                        \
         * Get an object derived from SetDataStorageBase to set values                             \
         */                                                                                        \
        return classtype<vartype>(&_##varname);                                                    \
    }

/*!
 * \brief Declares a named SetDataStorage getter function which can be controlled by a containing object
 * \param getname: The name of the getter function
 * \param varname: The name of the variable
 * \param classtype: The type of the class to be defined
 * \param vartype: The ctype of the variable
 * \param controller: A pointer to the controlling object
 */
#define TARDIGRADE_HYDRA_DECLARE_CONTROLLED_NAMED_SETDATASTORAGE_GETTER(getname, varname, classtype, vartype, \
                                                                        controller)                           \
    classtype<vartype> getname() {                                                                            \
        /*!                                                                                                   \
         * Get an object derived from SetDataStorageBase to set values                                        \
         */                                                                                                   \
        return classtype<vartype>(&_##varname, controller);                                                   \
    }

/*!
 * \brief Declares a SetDataStorage getter function
 * \param varname: The name of the variable
 * \param classtype: The type of the class to be defined
 * \param vartype: The ctype of the variable
 */
#define TARDIGRADE_HYDRA_DECLARE_SETDATASTORAGE_GETTER(varname, classtype, vartype) \
    TARDIGRADE_HYDRA_DECLARE_NAMED_SETDATASTORAGE_GETTER(get_SetDataStorage_##varname, varname, classtype, vartype)

/*!
 * \brief Declares a controlled SetDataStorage getter function
 * \param varname: The name of the variable
 * \param classtype: The type of the class to be defined
 * \param vartype: The ctype of the variable
 * \param controller: A pointer to the controlling object
 */
#define TARDIGRADE_HYDRA_DECLARE_CONTROLLED_SETDATASTORAGE_GETTER(varname, classtype, vartype, controller)            \
    TARDIGRADE_HYDRA_DECLARE_CONTROLLED_NAMED_SETDATASTORAGE_GETTER(get_SetDataStorage_##varname, varname, classtype, \
                                                                    vartype, controller)

/*!
 * \brief Declares a named getter function
 * \param getname: The name of the getter function
 * \param varname: The name of the variable
 * \param vartype: The ctype of the variable
 * \param uncall:  The function that is called if the variable is undefined.
 *     This should set the variable or throw an error.
 */
#define TARDIGRADE_HYDRA_DECLARE_NAMED_GETTER(getname, varname, vartype, uncall) \
    const vartype *getname() {                                                   \
        /*!                                                                      \
         * Get the value of varname                                              \
         */                                                                      \
        if (!_##varname.first) {                                                 \
            TARDIGRADE_ERROR_TOOLS_CATCH(uncall())                               \
        }                                                                        \
        return &_##varname.second;                                               \
    }

/*!
 * \brief Declares a getter function
 * \param varname: The name of the variable
 * \param vartype: The ctype of the variable
 * \param uncall:  The function that is called if the variable is undefined.
 *     This should set the variable or throw an error
 */
#define TARDIGRADE_HYDRA_DECLARE_GETTER(varname, vartype, uncall) \
    TARDIGRADE_HYDRA_DECLARE_NAMED_GETTER(get_##varname, varname, vartype, uncall)

/*!
 * \brief Declares a named setter function
 * \param setname: The name of the setter function
 * \param varname: The name of the variable
 * \param vartype: The ctype of the variable
 * \param setfun: The class member function that sets the variable
 */
#define TARDIGRADE_HYDRA_DECLARE_NAMED_SETTER(setname, varname, vartype, setfun) \
    void setname(const vartype &varname) {                                       \
        /*!                                                                      \
         * Set the storage variable varname of type vartype                      \
         *                                                                       \
         * \param &varname: The value of varname                                 \
         */                                                                      \
        TARDIGRADE_ERROR_TOOLS_CATCH(setfun(varname, _##varname))                \
    }

/*!
 * \brief Declares a setter function
 * \param varname: The name of the variable
 * \param vartype: The ctype of the variable
 * \param setfun: The class member function that sets the variable
 */
#define TARDIGRADE_HYDRA_DECLARE_SETTER(varname, vartype, setfun) \
    TARDIGRADE_HYDRA_DECLARE_NAMED_SETTER(set_##varname, varname, vartype, setfun)

/*!
 * \brief Declare a DataStorage class and the associated named setters and getters
 * \param context: The context (public, protected, private) that the macro should return to
 * \param setname: The name of the setter function
 * \param getname: The name of the getter function
 * \param getsetdatastoragename: The name of the getter function for the SetDataStorage object
 * \param varname: The name of the variable
 * \param classtype: The type of the class for the SetDataStorage object
 * \param vartype: The ctype of the variable
 * \param setfun: The class member function that sets the variable
 * \param uncall: The function that is called if the variable is not set
 */
#define TARDIGRADE_HYDRA_DECLARE_NAMED_STORAGE(context, setname, getname, getsetdatastoragename, varname, classtype, \
                                               vartype, setfun, uncall)                                              \
   private:                                                                                                          \
    DataStorage<vartype> _##varname;                                                                                 \
                                                                                                                     \
   public:                                                                                                           \
    TARDIGRADE_HYDRA_DECLARE_NAMED_SETTER(setname, varname, vartype, setfun)                                         \
   public:                                                                                                           \
    TARDIGRADE_HYDRA_DECLARE_NAMED_GETTER(getname, varname, vartype, uncall)                                         \
   public:                                                                                                           \
    TARDIGRADE_HYDRA_DECLARE_NAMED_SETDATASTORAGE_GETTER(getsetdatastoragename, varname, classtype, vartype)         \
    context:

/*!
 * \brief Declare a controlled DataStorage class and the associated named setters and getters
 * \param context: The context (public, protected, private) that the macro should return to
 * \param setname: The name of the setter function
 * \param getname: The name of the getter function
 * \param getsetdatastoragename: The name of the getter function for the SetDataStorage object
 * \param varname: The name of the variable
 * \param classtype: The type of the class for the SetDataStorage object
 * \param vartype: The ctype of the variable
 * \param setfun: The class member function that sets the variable
 * \param uncall: The function that is called if the variable is not set
 * \param controller: A pointer to the controlling object
 */
#define TARDIGRADE_HYDRA_DECLARE_CONTROLLED_NAMED_STORAGE(context, setname, getname, getsetdatastoragename, varname, \
                                                          classtype, vartype, setfun, uncall, controller)            \
   private:                                                                                                          \
    DataStorage<vartype> _##varname;                                                                                 \
                                                                                                                     \
   public:                                                                                                           \
    TARDIGRADE_HYDRA_DECLARE_NAMED_SETTER(setname, varname, vartype, setfun)                                         \
   public:                                                                                                           \
    TARDIGRADE_HYDRA_DECLARE_NAMED_GETTER(getname, varname, vartype, uncall)                                         \
   public:                                                                                                           \
    TARDIGRADE_HYDRA_DECLARE_CONTROLLED_NAMED_SETDATASTORAGE_GETTER(getsetdatastoragename, varname, classtype,       \
                                                                    vartype, controller)                             \
    context:

/*!
 * \brief Declare a DataStorage class and the associated setters and getters
 * \param context: The context (public, protected, private) that the macro should return to
 * \param varname: The name of the variable
 * \param vartype: The ctype of the variable
 * \param setfun: The class member function that sets the variable
 * \param uncall: The function that is called if the variable is not set
 */
#define TARDIGRADE_HYDRA_DECLARE_STORAGE(context, varname, vartype, setfun, uncall)                             \
    TARDIGRADE_HYDRA_DECLARE_NAMED_STORAGE(context, set_##varname, get_##varname, get_SetDataStorage_##varname, \
                                           varname, tardigradeHydra::SetDataStorageBase, vartype, setfun, uncall)

/*!
 * \brief Declare a named DataStorage class that uses setIterationData as the setter function
 * \param context: The context (public, protected, private) that the macro should return to
 * \param setname: The name of the setter
 * \param getname: The name of the getter function
 * \param varname: The name of the variable
 * \param vartype: The type of the variable
 * \param uncall: The function that is called if the variable is not set
 */
#define TARDIGRADE_HYDRA_DECLARE_NAMED_ITERATION_STORAGE(context, setname, getname, varname, vartype, uncall)      \
    TARDIGRADE_HYDRA_DECLARE_CONTROLLED_NAMED_STORAGE(context, setname, getname, get_SetDataStorage_##varname,     \
                                                      varname, SetDataStorageIteration, vartype, setIterationData, \
                                                      uncall, this)

/*!
 * \brief Declare a DataStorage class that uses setIterationData as the setter function
 * \param context: The context (public, protected, private) that the macro should return to
 * \param varname: The name of the variable
 * \param vartype: The type of the variable
 * \param uncall: The function that is called if the variable is not set
 */
#define TARDIGRADE_HYDRA_DECLARE_ITERATION_STORAGE(context, varname, vartype, uncall)                                 \
    TARDIGRADE_HYDRA_DECLARE_CONTROLLED_NAMED_STORAGE(context, set_##varname, get_##varname,                          \
                                                      get_SetDataStorage_##varname, varname, SetDataStorageIteration, \
                                                      vartype, setIterationData, uncall, this)

/*!
 * \brief Declare a DataStorage class that uses setNLStepData as the setter function
 * \param context: The context (public, protected, private) that the macro should return to
 * \param varname: The name of the variable
 * \param vartype: The type of the variable
 * \param uncall: The function that is called if the variable is not set
 */
#define TARDIGRADE_HYDRA_DECLARE_NLSTEP_STORAGE(context, varname, vartype, uncall)                                 \
    TARDIGRADE_HYDRA_DECLARE_CONTROLLED_NAMED_STORAGE(context, set_##varname, get_##varname,                       \
                                                      get_SetDataStorage_##varname, varname, SetDataStorageNLStep, \
                                                      vartype, setNLStepData, uncall, this)

/*!
 * \brief Declare a named DataStorage class that uses setPreviousData as the setter function
 * \param context: The context (public, protected, private) that the macro should return to
 * \param setname: The name of the setter
 * \param getname: The name of the getter function
 * \param varname: The name of the variable
 * \param vartype: The type of the variable
 * \param uncall: The function that is called if the variable is not set
 */
#define TARDIGRADE_HYDRA_DECLARE_NAMED_PREVIOUS_STORAGE(context, setname, getname, varname, vartype, uncall) \
    TARDIGRADE_HYDRA_DECLARE_NAMED_STORAGE(context, setname, getname, get_SetDataStorage_##varname, varname, \
                                           SetDataStoragePrevious, vartype, setPreviousData, uncall)

/*!
 * \brief Declare a DataStorage class that uses setPreviousData as the setter function
 * \param context: The context (public, protected, private) that the macro should return to
 * \param varname: The name of the variable
 * \param vartype: The type of the variable
 * \param uncall: The function that is called if the variable is not set
 */
#define TARDIGRADE_HYDRA_DECLARE_PREVIOUS_STORAGE(context, varname, vartype, uncall)                            \
    TARDIGRADE_HYDRA_DECLARE_NAMED_STORAGE(context, set_##varname, get_##varname, get_SetDataStorage_##varname, \
                                           varname, SetDataStoragePrevious, vartype, setPreviousData, uncall)

/*!
 * \brief Declare a named DataStorage class that uses setConstantData as the setter function
 * \param context: The context (public, protected, private) that the macro should return to
 * \param setname: The name of the setter
 * \param getname: The name of the getter function
 * \param varname: The name of the variable
 * \param vartype: The type of the variable
 * \param uncall: The function that is called if the variable is not set
 */
#define TARDIGRADE_HYDRA_DECLARE_NAMED_CONSTANT_STORAGE(context, setname, getname, varname, vartype, uncall) \
    TARDIGRADE_HYDRA_DECLARE_NAMED_STORAGE(context, setname, getname, get_setDataSTorage_##varname, varname, \
                                           SetDataStorageConstant, vartype, setConstantData, uncall)

/*!
 * \brief Declare a DataStorage class that uses setContantData as the setter function
 * \param context: The context (public, protected, private) that the macro should return to
 * \param varname: The name of the variable
 * \param vartype: The type of the variable
 * \param uncall: The function that is called if the variable is not set
 */
#define TARDIGRADE_HYDRA_DECLARE_CONSTANT_STORAGE(context, varname, vartype, uncall)                            \
    TARDIGRADE_HYDRA_DECLARE_NAMED_STORAGE(context, set_##varname, get_##varname, get_SetDataStorage_##varname, \
                                           varname, SetDataStorageConstant, vartype, setConstantData, uncall)

namespace tardigradeHydra {

    //    namespace core{

    /*!
     * A general class time for entries in a ValueStorage object
     */
    template <typename key_type, typename value_type>
    struct Entry {
        key_type   key;
        value_type value;
    };

    /*!
     * A compile-time map which co-locates related objects
     *
     * configuration is a class which defines the access to the object
     *
     * it must have a static member entries which is a static consexpr array of type
     * entry that defines key-value pairs as well as entries for key_type
     * and value_type which define the types of the keys and values
     * respectively
     */
    template <class configuration>
    class CompileTimeMap {
       public:
        /*!
         * Check if a key is in the map
         *
         * \param key: The key value to search in the map
         */
        constexpr bool key_in_map(typename configuration::key_type key) {
            for (const auto &entry : configuration::entries) {
                if (entry.key == key) {
                    return true;
                }
            }
            return false;
        }
        /*!
         * Get an entry to the map
         *
         * \param key: The key value to search in the map
         */
        // NOLINTBEGIN(clang-diagnostic-return-type)
        constexpr typename configuration::value_type get_value(typename configuration::key_type key) {
            for (const auto &entry : configuration::entries) {
                if (entry.key == key) {
                    return entry.value;
                }
            }
            // TODO: Figure out how to throw a compile time error
        }
        // NOLINTEND(clang-diagnostic-return-type)
    };

    /*!
     * Base class for data objects which defines the clear command
     */
    class dataBase {
       public:
        /*!
         * The function to erase the current values stored
         */
        virtual void clear() { TARDIGRADE_ERROR_TOOLS_CATCH(throw std::runtime_error("clear not implemented!")); }
    };

    /*!
     * Custom data storage object that allows for smart storage of objects
     */
    template <typename T>
    class DataStorage : public dataBase {
       public:
        bool first = false;  //!< The flag for whether the data has been stored

        T second = T();  //!< The stored data

        DataStorage() {};

        /*!
         * Constructor for a data-storage object setting first and second
         *
         * \param &_first: The flag for whether the data storage has been set
         * \param &_second: The data contained by the object
         */
        DataStorage(const bool &_first, const T &_second) : first(_first), second(_second) {}

        /*!
         * The function to erase the current values stored by setting first to false and clearing second
         */
        virtual void clear() {
            first = false;

            second.clear();
        }

        /*!
         * The function to set the current values to zero
         */
        virtual void zero() { std::fill(std::begin(second), std::end(second), 0); }

        /*!
         * The function to resize second and set the current values to zero
         *
         * \param size: The size to resize the vector to
         */
        virtual void zero(const unsigned int size) {
            second.resize(size);

            zero();
        }
    };

    /*!
     * The function to erase the current values stored by setting first to false and second to zero
     */
    template <>
    inline void DataStorage<int>::clear() {
        first = false;

        second = 0;
    }

    /*!
     * The function to set the value to zero
     */
    template <>
    inline void DataStorage<int>::zero() {
        second = 0;
    }

    /*!
     * The function to set the value to zero
     */
    template <>
    inline void DataStorage<int>::zero(const unsigned int size) {
        throw std::runtime_error("A scalar value cannot have a size!");
    }

    /*!
     * The function to erase the current values stored by setting first to false and second to zero
     */
    template <>
    inline void DataStorage<unsigned int>::clear() {
        first = false;

        second = 0;
    }

    /*!
     * The function to set the value to zero
     */
    template <>
    inline void DataStorage<unsigned int>::zero() {
        second = 0;
    }

    /*!
     * The function to set the value to zero
     */
    template <>
    inline void DataStorage<unsigned int>::zero(const unsigned int size) {
        throw std::runtime_error("A scalar value cannot have a size!");
    }

    /*!
     * A custom object that handles setting DataStorage objects in place
     * When the destructor is run the DataStorage is assumed to be set.
     */
    template <typename T>
    class SetDataStorageBase {
       public:
        //! Default constructor which points to no data
        SetDataStorageBase() : _ds(NULL) {}

        //! Constructor that reads in a data storage object and sets the value
        SetDataStorageBase(DataStorage<T> *ds) : _ds(ds) {
            if (_ds) {
                value = &_ds->second;
            }
        }

        //! Destructor which indicates that the data storage object has been set
        ~SetDataStorageBase() {
            if (_ds) {
                _ds->first = true;
            }
        }

        T *value;  //!< A pointer to the value of the data storage object

        //! Get the starting iterator of the data storage object
        auto begin() { return std::begin(*value); }

        //! Get the stopping iterator of the data storage object
        auto end() { return std::end(*value); }

        //! Initialize the data storage object to zero
        void zero() { _ds->zero(); }

        /*!
         * Resize (if possible) the data storage object and set to zero
         *
         * \param size: The size of the data storage object
         */
        void zero(const unsigned int size) { _ds->zero(size); }

        /*!
         * Create a zeroed matrix map of the quantity
         */
        template <typename X, unsigned int R, unsigned int C>
        auto zeroMap() {
            zero(R * C);
            return tardigradeHydra::getFixedSizeMatrixMap<X, R, C>(value->data());
        }

        /*!
         * Create a zeroed matrix map of the quantity with dynamic columns
         *
         * \param C: The number of columns in the matrix
         */
        template <typename X, unsigned int R>
        auto zeroMap(const unsigned int C) {
            zero(R * C);
            return tardigradeHydra::getDynamicColumnSizeMatrixMap<X, R>(value->data(), C);
        }

        /*!
         * Create a zeroed matrix map of the quantity as a dynamic-size matrix
         *
         * \param R: The number of rows in the matrix
         * \param C: The number of columns in the matrix
         */
        template <typename X>
        auto zeroMap(const unsigned int R, const unsigned int C) {
            zero(R * C);
            return tardigradeHydra::getDynamicSizeMatrixMap(value->data(), R, C);
        }

        /*!
         * Create a zeroed vector map of the quantity as a fixed size vector
         */
        template <typename X, unsigned int R>
        auto zeroMap() {
            zero(R);
            return tardigradeHydra::getFixedSizeVectorMap<X, R>(value->data());
        }

        /*!
         * Create a zeroed vector map of the quantity as a dynamic sized vector
         *
         * \param R: The number of rows in the vector
         */
        template <typename X>
        auto zeroMap(const unsigned int R) {
            zero(R);
            return tardigradeHydra::getDynamicSizeVectorMap(value->data(), R);
        }

       protected:
        DataStorage<T> *_ds;  //!< Pointer to the data for the data storage object
    };

    /*!
     * The base class for cached data which is to be reset at each iteration
     */
    template <typename container, typename T>
    class SetDataStorageIterationBase : public SetDataStorageBase<T> {
       public:
        /*!
         * Create a data storage object that will be reset at each new iteration
         *
         * \param *ds: The data storage object
         * \param *rp: The class that contains the data storage object. Must have a function addIterationData
         *     wich accepts a pointer to the data storage object
         */
        SetDataStorageIterationBase(DataStorage<T> *ds, container *rp) : SetDataStorageBase<T>(ds), _rp(rp) {}

        //! Destructor object that adds the data storage object to the iteration data list
        ~SetDataStorageIterationBase() { _rp->addIterationData(this->_ds); }

       protected:
        container *_rp;  //!< The containing hydraBase class
    };

    /*!
     * The base class for cached data which is to be reset at each non-linear step
     */
    template <typename container, typename T>
    class SetDataStorageNLStepBase : public SetDataStorageBase<T> {
       public:
        /*!
         * Create a data storage object that will be reset at each new nonlinear step
         *
         * \param *ds: The data storage object
         * \param *rp: The class that contains the data storage object. Must have a function addNLStepData
         *     wich accepts a pointer to the data storage object
         */
        SetDataStorageNLStepBase(DataStorage<T> *ds, container *rp) : SetDataStorageBase<T>(ds), _rp(rp) {}

        //! Destructor object that adds the data storage object to the nonlinear step data list
        ~SetDataStorageNLStepBase() { _rp->addNLStepData(this->_ds); }

       protected:
        container *_rp;  //!< The containing hydraBase class
    };

    // CoreDefinitions DataStorage specifications

    /*!
     * The function to erase the current values stored by setting first to false and second to zero
     */
    template <>
    inline void DataStorage<floatType>::clear() {
        first = false;

        second = 0;
    }

    /*!
     * The function to set the value to zero
     */
    template <>
    inline void DataStorage<floatType>::zero() {
        second = 0;
    }

    /*!
     * The function to set the value to zero
     */
    template <>
    inline void DataStorage<floatType>::zero(const unsigned int size) {
        throw std::runtime_error("A scalar value cannot have a size!");
    }

    //! A base class which defines functions for caching data
    class CachingDataBase {
       public:
        //! Constructor
        CachingDataBase() {}

        /*!
         * The function to add iteration data
         * All inheriting classes must implement this
         *
         * \param *data: The incoming data
         */
        virtual void addIterationData(dataBase *data) {
            TARDIGRADE_ERROR_TOOLS_CATCH(throw std::logic_error("addIterationData is not implemented"));
        }

        /*!
         * The function to add nonlinear step data
         * All inheriting classes must implement this
         *
         * \param *data: The incoming data
         */
        virtual void addNLStepData(dataBase *data) {
            TARDIGRADE_ERROR_TOOLS_CATCH(throw std::logic_error("addNLStepData is not implemented"));
        }

        /*!
         * Template function for adding iteration data. These values are cleared
         * every time the unknown vector is updated.
         *
         * \param &data: The data to be added
         * \param &storage: The storage to add the data to
         */
        template <class T>
        void setIterationData(const T &data, DataStorage<T> &storage) {
            storage.second = data;

            storage.first = true;

            addIterationData(&storage);
        }

        /*!
         * Template function for adding nonlinear step data. These values are cleared
         * every time the nonlinear step is advanced.
         *
         * \param &data: The data to be added
         * \param &storage: The storage to add the data to
         */
        template <class T>
        void setNLStepData(const T &data, DataStorage<T> &storage) {
            storage.second = data;

            storage.first = true;

            addNLStepData(&storage);
        }

        /*!
         * Template function for adding previous data
         *
         * \param &data: The data to be added
         * \param &storage: The storage to add the data to
         */
        template <class T>
        void setPreviousData(const T &data, DataStorage<T> &storage) {
            storage.second = data;

            storage.first = true;
        }

        /*!
         * Template function for adding constant data
         *
         * \param &data: The data to be added
         * \param &storage: The storage to add the data to
         */
        template <class T>
        void setConstantData(const T &data, DataStorage<T> &storage) {
            storage.second = data;

            storage.first = true;
        }

       protected:
        /*!
         * Class which defines data storage objects which are reset at each iteration
         */
        template <typename T>
        class SetDataStorageIteration : public SetDataStorageIterationBase<CachingDataBase, T> {
           public:
            using tardigradeHydra::SetDataStorageIterationBase<CachingDataBase, T>::SetDataStorageIterationBase;
        };

        /*!
         * Class which defines data storage objects which are reset at each nonlinear step
         */
        template <typename T>
        class SetDataStorageNLStep : public SetDataStorageNLStepBase<CachingDataBase, T> {
           public:
            using tardigradeHydra::SetDataStorageNLStepBase<CachingDataBase, T>::SetDataStorageNLStepBase;
        };

        //! Class which defines data storage objects for values defined at the previous timestep
        template <typename T>
        class SetDataStoragePrevious : public SetDataStorageBase<T> {
           public:
            /*!
             * Constructor for data storage objects for temporally previous objects
             *
             * \param *ds: The data storage object to modify
             */
            SetDataStoragePrevious(DataStorage<T> *ds) : SetDataStorageBase<T>(ds) {}
        };

        /*!
         * Class that is a constant data storage object
         */
        template <typename T>
        class SetDataStorageConstant : public SetDataStorageBase<T> {
           public:
            /*!
             * Constructor for constant data storage objects
             *
             * \param *ds: The data storage object
             */
            SetDataStorageConstant(DataStorage<T> *ds) : SetDataStorageBase<T>(ds) {}
        };
    };

    //    }

}  // namespace tardigradeHydra

#include "tardigrade_SetDataStorage.cpp"
#include "tardigrade_SetDataStorage.tpp"

#endif
