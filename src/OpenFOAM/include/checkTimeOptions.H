// Deprecated include (2024-05) - prefer timeSelector
// ----------------------------------------------------------------------------

Foam::label startTime = 0;

// Unless -constant is present, skip startTime if it is "constant"
#include "checkConstantOption.H"

// Check -time and -latestTime options
#include "checkTimeOption.H"

// ************************************************************************* //
