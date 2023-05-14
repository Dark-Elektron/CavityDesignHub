Set WshShell = CreateObject("WScript.Shell")
Set fso = CreateObject("Scripting.FileSystemObject")
WScript.Echo "hello, world"

' Get current path
sCurPath = fso.GetAbsolutePathName(".")
sParentPath = fso.GetParentFolderName(WScript.ScriptFullName)

WScript.Echo sCurPath
WScript.Echo sParentPath

'Get the command line arguments passed by DAKOTA, so we work on the correct files
paramsfile  = WScript.Arguments.Item(0)
resultsfile = WScript.Arguments.Item(1)

' Run the Python script and capture the output
If (fso.FileExists(sParentPath & "\cubature_nodes.xlsx")) Then
	cmd = "python " & sParentPath & "\py_dakota.py " & sCurPath & "\" & paramsfile & " " & sCurPath & "\" & resultsfile & " 2" & " ALL"
Else
	cmd = "python " & sParentPath & "\py_dakota.py " & sCurPath & "\" & paramsfile & " " & sCurPath & "\" & resultsfile & " 2" & " ONLY_NODES"
End If

WScript.Echo cmd

Set exec = WshShell.Exec(cmd)
output = exec.StdOut.ReadAll()

' Print the output in the command prompt
WScript.Echo "Python Output: " & output
