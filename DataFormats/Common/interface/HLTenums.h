#ifndef DataFormats_Common_HLTenums_h
#define DataFormats_Common_HLTenums_h

/** \brief HLT enums
 *
 *  Definition of common HLT enums
 *
 *  $Date: 2013/01/24 09:48:36 $
 *  $Revision: 1.1 $
 *
 *  \author Martin Grunewald
 *
 */

namespace edm
{
  namespace hlt
    {

      /// status of a trigger path
      enum HLTState {Ready=0,     ///< not [yet] run
		     Pass =1,     ///< accept
		     Fail =2,     ///< reject
		     Exception=3  ///< error
      };

    }
}

#endif // DataFormats_Common_HLTenums_h
