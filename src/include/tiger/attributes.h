#ifndef _ATTRIBUTES_H_
#define _ATTRIBUTES_H_

#pragma once

#include <unordered_map>
#include <vector>

const std::string ATTRIB_BOTTOM_POINT = "bottom_point";
const std::string ATTRIB_BOUNDARY = "boundary";
const std::string ATTRIB_CENTER = "center";
const std::string ATTRIB_CENTROID = "centroid";
const std::string ATTRIB_COLOR = "color";
const std::string ATTRIB_COPLANAR = "coplanar";
const std::string ATTRIB_DIRECTION_VALUE = "dir_value";
const std::string ATTRIB_DIRECTION_VECTOR = "dir_vector";
const std::string ATTRIB_DUAL_INDEX = "dual_index";
const std::string ATTRIB_ENQUEUED = "enqueued";
const std::string ATTRIB_EQUILIBRIUM_INDEX = "equillibrium_index";
//const std::string ATTRIB_FACE_INDEX = "face_index";
const std::string ATTRIB_FIND_AT = "find_at";
const std::string ATTRIB_FIXED = "fixed";
const std::string ATTRIB_HALFEDGE = "halfedge";
const std::string ATTRIB_ID = "id";
const std::string ATTRIB_INDEX = "index";
const std::string ATTRIB_IN_EQUILIBRIUM = "in_equilibrum";
const std::string ATTRIB_LOCATION = "location";
const std::string ATTRIB_MIDPOINT = "midpoint";
const std::string ATTRIB_NEGATIVE_INDEX = "negative_index";
const std::string ATTRIB_NEGATIVE_BLOCK_INDEX = "negative_block_index";
const std::string ATTRIB_NORMAL = "normal";
const std::string ATTRIB_ORIGINAL = "original";
//const std::string ATTRIB_BLOCK_INDEX = "block_index";
//const std::string ATTRIB_BLOCK_ENABLED = "block_enabled";
const std::string ATTRIB_POSITIVE_INDEX = "positive_index";
const std::string ATTRIB_POSITIVE_BLOCK_INDEX = "positive_block_index";
const std::string ATTRIB_ROTATED_VECTOR = "rotated_vector";
const std::string ATTRIB_SPLIT = "split";
const std::string ATTRIB_TOP_POINT = "top_point";
const std::string ATTRIB_MODIFIED = "modified";
const std::string ATTRIB_VISITED = "visited";

/*
The base class for a dynamic attribute. It literally has nothing, but we require it.
*/
class Attr
{

public:

    /*
    Destructor of the class.
    */
    virtual ~Attr() {}

};


/*
The class representing a dynamic attribute.
*/
template <typename T>
class Attribute : public Attr
{

private:

    // The data of the attribute
    T m_data;

public:

    /*
    Constructor of the class.
    @param T d The data of the attribute.
    */
    Attribute(T d) : m_data(d) {}

    /*
    Destructor of the class
    */
    virtual ~Attribute() {}

    /*
    Get the data of the attribute.
    @return T The data of the attribute.
    */
    T data()
    {
        return m_data;
    }

    /*
    Set the data of the attribute.
    @param T d The data for the attribute.
    */
    void data(T d)
    {
        m_data = d;
    }

};


/*
The class representing a set of dynamic attributes.
*/
class DynamicAttributes
{

private:

    // The map that stores the attributes
    std::unordered_map<std::string, Attr *> m_attrs;

public:

    /*
    Constructor of the class.
    */
    DynamicAttributes()
    {
    }

    /*
    Destructor of the class.
    */
    ~DynamicAttributes()
    {
        Clear();
    }

    /*
    Clears the dynamic attributes object. All of them are deleted.
    */
    void Clear()
    {
        // Get the keys of the stored attributes
        std::vector<std::string> keys;
        Names(keys);

        // Traverse through the keys and delete their content
        for (auto it = keys.begin(); it != keys.end(); ++it) 
        {
            Erase(*it);
        }

        // Clear the attributes map
        m_attrs.clear();
    }

    /*
    Erases the attribute from the map
    @param std::string name The attribute name
    */
    void Erase(const std::string & name)
    {
        // Search for the attribute name in the attributes map
        auto it = m_attrs.find(name);

        // If the attribute name was found then erase it from the map
        if (it != m_attrs.end())
        {
            // Delete the object containing the value. Then, delete the entry from the map
            delete it->second;
            m_attrs.erase(it);
        }
    }

    /*
    Get the value of the indicated attribute.
    @param T & data The reference to the attribute value.
    @return bool Indicates if the attribute was found.
    */
    template <typename T>
    bool Get(const std::string & name, T & data) const 
    {
        // Find the attribute in the map
        auto it = m_attrs.find(name);

        // If the attribute is in the map then get its value and return true; otherwise, return 
        // false
        if (it != m_attrs.end())
        {
            Attribute<T> * p = (Attribute<T> *)m_attrs.at(name);
            data = p->data();
            return true;
        }
        else
        {
            return false;
        }
    }

    /*
    Checks whether the attribute name is registered.
    @param const std::string name The name of the attribute.
    @return bool true if the name is found, otherwise false.
    */
    bool Has(const std::string & name) const 
    {
        // If the attribute name was found then return true, otherwise return false
        return (m_attrs.find(name) != m_attrs.end());
    }

    /*
    Populates a vector with the names of the stored attributes.
    @param std::vector<std::string> & names The reference to the vector for the attribute names.
    @return bool Indicates if the vector was populated.
    */
    bool Names(std::vector<std::string> & names) const 
    {
        size_t nKeys = m_attrs.size();

        names.clear();
        names.resize(nKeys);

        size_t i = 0;

        // Traverse through the attributes map and push the keys into the names vector
        for (auto it = m_attrs.begin(); it != m_attrs.end(); ++it)
        {
            names[i++] = it->first;
        }

        return true;
    }

    /*
    @param const std::string & name
    @param T data
    */
    template <typename T>
    void Set(const std::string & name, T data)
    {
        // Search for the attribute name in the attributes map
        auto it = m_attrs.find(name);

        // If the name was not found then add the new attribute to the map; otherwise, delete the 
        // current attribute and make a new one with the given data
        if (it == m_attrs.end())
        {
            m_attrs.insert(std::pair<std::string, Attr *>(name, new Attribute<T>(data)));
        }
        else
        {
            delete it->second;
            it->second = new Attribute<T>(data);
        }
    }

    /*
    */
    template <typename T>
    T operator [](std::string & key) const
    {
        // Search for the attribute name in the attributes map
        auto it = m_attrs.find(key);
        assert(key != m_attrs.end());

        return it->second;
    }

    /*
    */
    template <typename T>
    T & operator [](std::string & key)
    {
        // Search for the attribute name in the attributes map
        auto it = m_attrs.find(key);

        // If the name was not found then add the new attribute to the map; otherwise, delete the 
        // current attribute and make a new one with the given data
        if (it == m_attrs.end())
        {
            m_attrs.insert(std::pair<std::string, Attr *>(key, new Attribute<T>(data)));
        }
        else
        {
            delete it->second;
            it->second = new Attribute<T>(data);
        }

        return it->second;
    }

};

/*
*/
class AttributedClass 
{

protected:

    // The dynamic attributes
    DynamicAttributes m_attributes;

public:

    /*
    Constructor of the class.
    */
    AttributedClass() {};

    /*
    Destructor of the class.
    */
    ~AttributedClass()
    {
        m_attributes.Clear();
    };

    /*
    Returns the reference to the dynamic attributes.
    @return Attributes & The reference to the dynamic attributes.
    */
    DynamicAttributes & Attributes() 
    {
        return m_attributes;
    };

};

#endif
