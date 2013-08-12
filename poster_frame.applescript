on run (arguments)
	set pwd to do shell script "pwd"
	log (do shell script ("rm -f " & (item 2 of arguments)))
	set inFile to POSIX file (pwd & "/" & item 1 of arguments) as text
	set outFile to POSIX file (pwd & "/" & item 2 of arguments) as text
	log inFile
	log outFile
	tell application "QuickTime Player 7"
		set d to open file inFile
		select none d
		copy d
		close d
		set d2 to make new document
		paste d2
		export d2 to outFile as image sequence using settings "BMP, 25 fps"
		close d2
	end tell
end run
