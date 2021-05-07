#ifndef _MESSAGE_H_
#define _MESSAGE_H_

#pragma once

//
// The class representing a message. Instances of this class are usually used for indicating missing 
// requirements for the procedures in a workspace.
//
class Message 
{

public:

    // The different message types
    enum MESSAGE_TYPE
    {
        MSG_UNKNOWN,
        MSG_BLOCKS_MODIFIED,
        MSG_NO_BLOCKS,
        MSG_NO_EDGE_DIRECTIONS,
        MSG_NO_EVEN_SIDED_TILES,
        MSG_NO_INTERFACES,
        MSG_NO_TESSELLATION,
        MSG_NO_TILE_CENTERS,
        MSG_OK,
        MSG_DEFAULT
    };

private:

    // 
    const char * STR_UNKNOWN                = "Unknown message type.";
    const char * STR_BLOCKS_MODIFIED        = "blocks are modified.";
    const char * STR_NO_BLOCKS              = "There are no blocks.";
    const char * STR_NO_EDGE_DIRECTIONS     = "Edges don't have directions.";
    const char * STR_NO_EVEN_SIDED_TILES    = "Tiles don't have an even number of sides.";
    const char * STR_NO_INTERFACES          = "There are no interface polygons.";
    const char * STR_NO_TESSELLATION        = "There is no tessellation";
    const char * STR_NO_TILE_CENTERS        = "Tiles don't have centers.";
    const char * STR_OK                     = "OK";
    const char * STR_DEFAULT                = "Default message";

    // The type of the message
    MESSAGE_TYPE m_type;

public:

	//
    // Constructor of the class. Object initializes with MSG_DEFAULT type.
	//
	Message();

    //
    // Constructor of the class.
    // @param MESSAGE_TYPE type The type of the message.
    //
    Message(MESSAGE_TYPE type);

    //
    // Destructor of the class.
    //
    ~Message();

    //
    // Returns the message text.
    // @return const char * The text of the message.
    //
    const char * getText() const;

    //
    // Returns the message type.
    // @return MESSAGE_TYPE The type of the message.
    //
    MESSAGE_TYPE getType() const;

    //
    // Sets the content of the message based on the given message type.
    // @param MESSAGE_TYPE type The type of the message.
    //
    void set(MESSAGE_TYPE type);

};

#endif
