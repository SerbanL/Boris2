#pragma once

#include "CompileFlags.h"
#if GRAPHICS == 1

#include <string>

#include "BorisLib.h"

//Simulation object holds state and action handlers
class TextObject;
class Simulation;
typedef Simulation HandlerObject;

//Interactive object states enum

//IOS_DELETING : State handler has marked this TextObject as needing deletion : caller must implement this
//IOS_DELETINGPARAGRAPH : State handler has marked the entire paragraph containing the TextObject as needing deletion : caller must implement this
//IOS_REPLACINGPARAGRAPH : State handler has marked the entire paragraph containing the TextObject as needing replaceing : caller must implement this
//
//IOS_ISLASTINLIST, IOS_WASLASTINLIST : used to handle representations of a list of things through action objects (e.g. a list of simulation stages), allowing the list to be updated automatically when elements are added at the end of the list.
//IOS_MISALIGNEDTAO : A text alignment object (tao) has detected misalignment between the ending rectangle coordinates of synchronised taos (of which it is one). They all need to be adjusted at the same time (spaces removed or added in their text)

enum IOS_ {IOS_ON = 0, IOS_OFF, IOS_ISLASTINLIST, IOS_WASLASTINLIST, IOS_DELETING, IOS_DELETINGPARAGRAPH, IOS_REPLACINGPARAGRAPH, IOS_MISALIGNEDTAO};

//////////////////////////////////////////////////////////////////////////////////////////////////
//Wrapper for action outcome codes (AO_) with any returned text - objects of this type are returned by the action handler

struct InteractiveObjectActionOutcome {

	//the action outcome code
	int actionOutcome;

	//any returned text
	std::string text;

	InteractiveObjectActionOutcome(void) : actionOutcome(0) {}
	InteractiveObjectActionOutcome(int actionOutcome_) : actionOutcome(actionOutcome_) {}

	//assignment operator
	InteractiveObjectActionOutcome& operator=(int actionOutcome_) { actionOutcome = actionOutcome_; return *this; }

	bool operator==(int actionOutcome_) const { return actionOutcome == actionOutcome_; }
	bool operator!=(int actionOutcome_) const { return actionOutcome != actionOutcome_; }
};

//////////////////////////////////////////////////////////////////////////////////////////////////
// TextObject -> Simulation Functionoid

//Wrapper for Simulation function member pointers, specifically used for interactive TextObjects. To use it, define a SimTOFunct pointer and initialise it in Simulation constructor.
//Two methods needed: one which defines an action handler (ActionHandler) and another which defines a state handler (StateHandler).
//ActionHandler implements functionality when the TextObject is interacted with. Response depends on TextObject type (callerId), action code (e.g. left mouse click), its state and text identifier.
//ActionHandler does not change any object properties, it simply implements functionality depending on stored object properties at the time it was interacted with.
//
//StateHandler can change state, textId, and the way the TextObject looks depending on its state - this is needed since the TextObject state can change due to actions other than interacting with it, thus StateHandler is called before the object is drawn.
//Return true if a change was made. This lets the caller know text lines must be recalculated (StateHandler is drawn by the top Draw method in TextLines).
//Special IOS_ states must be handled by the caller also (see IOS_ enum description)

//InteractiveObjectStateChange returned by state handler method : pass message to caller
struct InteractiveObjectStateChange {

	//flag to signal this text object has changed state : check its IOS_ state and take appropriate action.
	bool stateChanged;

	//some IOS_ states rely on a text message returned by the state handler
	std::string textMessage;

	InteractiveObjectStateChange(void) : stateChanged(false) {}

	//assignment operator : set state change flag
	InteractiveObjectStateChange& operator=(bool stateChanged_) { stateChanged = stateChanged_; return *this; }

	InteractiveObjectStateChange& operator()(bool stateChanged_) { stateChanged = stateChanged_; return *this; }

	bool operator==(const bool& rhs) const { return stateChanged == rhs; }
	bool operator!=(const bool& rhs) const { return stateChanged != rhs; }
	bool operator!(void) const { return !stateChanged; }
};

struct InteractiveObjectProperties {

public:

	//text object type. type has major and minor entries. These are typically defined in an enum and are used to decide what action is taken when an object is interacted with
	int majorId, minorId;

	//store a value from the IOS_ enum: checked by action handler and modified if needed by state handler
	IOS_ state;

	//auxiliary number identifier
	int auxId;

	//a text identifier for this object
	std::string textId;

	//used to interact another object with this - set this Id before calling the Action handler with AC_INTERACTOBJECTS code
	INT2 interactingObjectId;

public:

	InteractiveObjectProperties(void) 
	{
		majorId = 0;
		minorId = 0;
		state = IOS_ON;
		auxId = 0;
		textId = "";
		interactingObjectId = INT2();
	}

	InteractiveObjectProperties(int majorId, int minorId = 0, int auxId = 0, std::string textId = "") 
	{
		this->majorId = majorId;
		this->minorId = minorId;
		state = IOS_ON;
		this->auxId = auxId;
		this->textId = textId;
		interactingObjectId = INT2();
	}
};

//Define function pointer to member function types so they can be used throughout the program. The actual classes are defined elsewhere.

enum AO_;

typedef InteractiveObjectActionOutcome (HandlerObject::*SimAH)(int actionCode, InteractiveObjectProperties& iop, TextObject *pTO);
typedef InteractiveObjectStateChange (HandlerObject::*SimSH)(InteractiveObjectProperties& iop, TextObject *pTO);

class SimTOFunct {

	//Action and state handler function pointers, together with the class object (Simulation class) they are held in
	HandlerObject *pHandler;

	SimAH ActionHandler;
	SimSH StateHandler;

public:

	//properties of the interactive object storing a SimTOFunct pointer.
	InteractiveObjectProperties iop;

public:
	
	SimTOFunct(HandlerObject *pHandler, SimAH ActionHandler, SimSH StateHandler) 
	{
		this->pHandler = pHandler;
		this->ActionHandler = ActionHandler;
		this->StateHandler = StateHandler;
	}
	
	//copy constructor when used with new operator
	SimTOFunct(const SimTOFunct *pSTF) 
	{
		pHandler = pSTF->pHandler;
		ActionHandler = pSTF->ActionHandler;
		StateHandler = pSTF->StateHandler;
		iop = pSTF->iop;
	}

	~SimTOFunct() {}

	void set_iop(InteractiveObjectProperties iop) { this->iop = iop; }

	//Calls to Action and State handlers
	InteractiveObjectActionOutcome Action(int actionCode, TextObject *pTO) { return CALLFP(pHandler, ActionHandler)(actionCode, iop, pTO); }
	InteractiveObjectStateChange State(TextObject *pTO) { return CALLFP(pHandler, StateHandler)(iop, pTO); }
};

#endif